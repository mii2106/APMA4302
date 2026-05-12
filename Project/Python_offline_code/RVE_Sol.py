from pathlib import Path

import gmsh
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem
from dolfinx.io import gmshio, XDMFFile
from dolfinx_mpc import (MultiPointConstraint,assemble_matrix,assemble_vector,apply_lifting,)


# -----------------------------------------------------------------------------
# MPI / output
# -----------------------------------------------------------------------------
world_comm = MPI.COMM_WORLD
rank = world_comm.rank
size = world_comm.size
t_total_start = MPI.Wtime()

# Task-parallel mode: every MPI rank owns a complete copy of the RVE mesh and
# solves a subset of the independent EBH right-hand sides. This is better than
# domain decomposition for small meshes with many partition problems.
comm = MPI.COMM_SELF
Path("roh_cache").mkdir(exist_ok=True)
Path("vtk_outputs").mkdir(exist_ok=True)


# -----------------------------------------------------------------------------
# Geometry / mesh parameters
# -----------------------------------------------------------------------------
Lx = Ly = Lz = 1.0
fv = 0.5
R = np.sqrt(fv * (Lx * Ly * Lz) / (np.pi * Lz))
cx, cy = Lx / 2.0, Ly / 2.0
h = 0.30

# Physical tags
TAG_MATRIX = 1
TAG_FIBER = 2


# -----------------------------------------------------------------------------
# Material / EBH user parameters
# -----------------------------------------------------------------------------
E_matrix = 100.0
nu_matrix = 0.3
E_fiber = 10000.0
nu_fiber = 0.3

PART_TYPE = "angular"
M_PARTS = 32
STRIP_AXIS = 3
FIBER_AXIS = 3
PART_FILE = "partitions.dat"
AREA_FILE = "area.dat"

PART_TYPE = "radial_angular"
N_ANGULAR = 24
N_RADIAL = 4
M_PARTS = N_ANGULAR * N_RADIAL
FIBER_AXIS = 3

NV = 6
VOIGT_LABELS = np.array(["xx", "yy", "zz", "xy", "yz", "xz"])


gmsh.initialize()
gmsh.open("modelo_rve.msh")
domain, cell_tags, facet_tags = gmshio.model_to_mesh(gmsh.model,comm,0,gdim=3)
gmsh.finalize()


# -----------------------------------------------------------------------------
# Function spaces
# -----------------------------------------------------------------------------
gdim = domain.geometry.dim
V = fem.functionspace(domain, ("Lagrange", 1, (gdim,)))
W = fem.functionspace(domain, ("Discontinuous Lagrange", 0))


# -----------------------------------------------------------------------------
# Material fields, piecewise constant per cell
# -----------------------------------------------------------------------------
fib_cells = cell_tags.find(TAG_FIBER)
mat_cells = cell_tags.find(TAG_MATRIX)

E_fun = fem.Function(W)
nu_fun = fem.Function(W)
E_fun.name = "E"
nu_fun.name = "nu"

E_fun.x.array[:] = E_matrix
nu_fun.x.array[:] = nu_matrix
E_fun.x.array[fib_cells] = E_fiber
nu_fun.x.array[fib_cells] = nu_fiber

mu_fun = fem.Function(W)
lmbda_fun = fem.Function(W)
mu_fun.name = "mu"
lmbda_fun.name = "lambda"

mu_fun.x.array[:] = E_fun.x.array / (2.0 * (1.0 + nu_fun.x.array))
lmbda_fun.x.array[:] = (E_fun.x.array * nu_fun.x.array) / (
    (1.0 + nu_fun.x.array) * (1.0 - 2.0 * nu_fun.x.array)
)

I = ufl.Identity(gdim)


def eps(u):
    return ufl.sym(ufl.grad(u))


def C_eps(e):
    return 2.0 * mu_fun * e + lmbda_fun * ufl.tr(e) * I


# -----------------------------------------------------------------------------
# Periodic boundary constraints using dolfinx_mpc
# This is the periodic block from your script, kept essentially unchanged.
# -----------------------------------------------------------------------------
mpc = MultiPointConstraint(V)
tol = 1.0e-5
near = lambda a, b: np.isclose(a, b, atol=tol)

xmin, xmax = -0.5, 0.5
ymin, ymax = -0.5, 0.5
zmin, zmax = 0.0, 1.0

Lx = xmax - xmin
Ly = ymax - ymin
Lz = zmax - zmin

corners_mask = lambda x: (
    (near(x[0], xmin) & near(x[1], ymin) & near(x[2], zmin))
    | (near(x[0], xmax) & near(x[1], ymin) & near(x[2], zmin))
    | (near(x[0], xmin) & near(x[1], ymax) & near(x[2], zmin))
    | (near(x[0], xmax) & near(x[1], ymax) & near(x[2], zmin))
    | (near(x[0], xmin) & near(x[1], ymin) & near(x[2], zmax))
    | (near(x[0], xmax) & near(x[1], ymin) & near(x[2], zmax))
    | (near(x[0], xmin) & near(x[1], ymax) & near(x[2], zmax))
    | (near(x[0], xmax) & near(x[1], ymax) & near(x[2], zmax))
)

pin_dofs = fem.locate_dofs_geometrical(V, corners_mask)
u_zero = fem.Function(V)
u_zero.x.array[:] = 0.0
bc_pin = fem.dirichletbc(u_zero, pin_dofs)
bcs = [bc_pin]

# Face constraints, excluding edges
ind_xmin_faces = lambda x: near(x[0], xmin) & ~(
    near(x[1], ymin) | near(x[1], ymax) | near(x[2], zmin) | near(x[2], zmax)
)
ind_ymin_faces = lambda x: near(x[1], ymin) & ~(
    near(x[0], xmin) | near(x[0], xmax) | near(x[2], zmin) | near(x[2], zmax)
)
ind_zmin_faces = lambda x: near(x[2], zmin) & ~(
    near(x[0], xmin) | near(x[0], xmax) | near(x[1], ymin) | near(x[1], ymax)
)

rel_x = lambda x: np.vstack((x[0] + Lx, x[1], x[2]))
rel_y = lambda x: np.vstack((x[0], x[1] + Ly, x[2]))
rel_z = lambda x: np.vstack((x[0], x[1], x[2] + Lz))

mpc.create_periodic_constraint_geometrical(V, ind_xmin_faces, rel_x, bcs, 1.0)
mpc.create_periodic_constraint_geometrical(V, ind_ymin_faces, rel_y, bcs, 1.0)
mpc.create_periodic_constraint_geometrical(V, ind_zmin_faces, rel_z, bcs, 1.0)

# Edges parallel to z
edges_xmin_ymin = lambda x: near(x[0], xmin) & near(x[1], ymin) & ~(
    near(x[2], zmin) | near(x[2], zmax)
)
rel_xmin_ymin = lambda x: np.vstack((x[0] + Lx, x[1] + Ly, x[2]))
mpc.create_periodic_constraint_geometrical(V, edges_xmin_ymin, rel_xmin_ymin, bcs, 1.0)

edges_xmin_ymax = lambda x: near(x[0], xmin) & near(x[1], ymax) & ~(
    near(x[2], zmin) | near(x[2], zmax)
)
rel_xmin_ymax = lambda x: np.vstack((x[0] + Lx, x[1], x[2]))
mpc.create_periodic_constraint_geometrical(V, edges_xmin_ymax, rel_xmin_ymax, bcs, 1.0)

edges_xmax_ymin = lambda x: near(x[0], xmax) & near(x[1], ymin) & ~(
    near(x[2], zmin) | near(x[2], zmax)
)
rel_xmax_ymin = lambda x: np.vstack((x[0], x[1] + Ly, x[2]))
mpc.create_periodic_constraint_geometrical(V, edges_xmax_ymin, rel_xmax_ymin, bcs, 1.0)

# Edges parallel to y
edges_xmin_zmin = lambda x: near(x[0], xmin) & near(x[2], zmin) & ~(
    near(x[1], ymin) | near(x[1], ymax)
)
rel_xmin_zmin = lambda x: np.vstack((x[0] + Lx, x[1], x[2] + Lz))
mpc.create_periodic_constraint_geometrical(V, edges_xmin_zmin, rel_xmin_zmin, bcs, 1.0)

edges_xmin_zmax = lambda x: near(x[0], xmin) & near(x[2], zmax) & ~(
    near(x[1], ymin) | near(x[1], ymax)
)
rel_xmin_zmax = lambda x: np.vstack((x[0] + Lx, x[1], x[2]))
mpc.create_periodic_constraint_geometrical(V, edges_xmin_zmax, rel_xmin_zmax, bcs, 1.0)

edges_xmax_zmin = lambda x: near(x[0], xmax) & near(x[2], zmin) & ~(
    near(x[1], ymin) | near(x[1], ymax)
)
rel_xmax_zmin = lambda x: np.vstack((x[0], x[1], x[2] + Lz))
mpc.create_periodic_constraint_geometrical(V, edges_xmax_zmin, rel_xmax_zmin, bcs, 1.0)

# Edges parallel to x
edges_ymin_zmin = lambda x: near(x[1], ymin) & near(x[2], zmin) & ~(
    near(x[0], xmin) | near(x[0], xmax)
)
rel_ymin_zmin = lambda x: np.vstack((x[0], x[1] + Ly, x[2] + Lz))
mpc.create_periodic_constraint_geometrical(V, edges_ymin_zmin, rel_ymin_zmin, bcs, 1.0)

edges_ymin_zmax = lambda x: near(x[1], ymin) & near(x[2], zmax) & ~(
    near(x[0], xmin) | near(x[0], xmax)
)
rel_ymin_zmax = lambda x: np.vstack((x[0], x[1] + Ly, x[2]))
mpc.create_periodic_constraint_geometrical(V, edges_ymin_zmax, rel_ymin_zmax, bcs, 1.0)

edges_ymax_zmin = lambda x: near(x[1], ymax) & near(x[2], zmin) & ~(
    near(x[0], xmin) | near(x[0], xmax)
)
rel_ymax_zmin = lambda x: np.vstack((x[0], x[1], x[2] + Lz))
mpc.create_periodic_constraint_geometrical(V, edges_ymax_zmin, rel_ymax_zmin, bcs, 1.0)

mpc.finalize()


# -----------------------------------------------------------------------------
# Voigt/tensor helpers
# Convention: voigt = [xx, yy, zz, xy, yz, xz]
# Off-diagonal entries in tensor are 0.5 for unit engineering shear.
# -----------------------------------------------------------------------------
def unit_E(l: int):
    E = np.zeros((3, 3), dtype=PETSc.ScalarType)
    if l == 0:
        E[0, 0] = 1.0
    elif l == 1:
        E[1, 1] = 1.0
    elif l == 2:
        E[2, 2] = 1.0
    elif l == 3:
        E[0, 1] = 0.5
        E[1, 0] = 0.5
    elif l == 4:
        E[1, 2] = 0.5
        E[2, 1] = 0.5
    elif l == 5:
        E[0, 2] = 0.5
        E[2, 0] = 0.5
    else:
        raise ValueError("l must be in [0, 5]")
    return E


def to_voigt(M):
    return np.array(
        [M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[1, 2], M[0, 2]], dtype=float
    )


def masked_tensor(E, chi):
    return ufl.as_tensor(
        [
            [chi * E[0, 0], chi * E[0, 1], chi * E[0, 2]],
            [chi * E[1, 0], chi * E[1, 1], chi * E[1, 2]],
            [chi * E[2, 0], chi * E[2, 1], chi * E[2, 2]],
        ]
    )


def average_strain_over_mask(u_sol, chi, volume, add_macro=None):
    """Return <eps(u_sol) + add_macro> over one partition mask."""
    strain = eps(u_sol)
    if add_macro is not None:
        strain = strain + ufl.as_tensor(add_macro)

    M = np.zeros((3, 3), dtype=float)
    dx = ufl.dx(domain=domain)
    for i in range(3):
        for j in range(3):
            local_val = fem.assemble_scalar(fem.form(chi * strain[i, j] * dx))
            M[i, j] = comm.allreduce(local_val, op=MPI.SUM) / volume
    return to_voigt(M)


# -----------------------------------------------------------------------------
# EBH partitions
# Follows ebh_offline_3d.f90:
#   part_id = 1      -> all fiber cells
#   part_id = 2..M+1 -> M matrix strips/sectors
# -----------------------------------------------------------------------------
def _axis_index(axis: int) -> int:
    if axis not in (1, 2, 3):
        raise ValueError("axis must be 1, 2, or 3")
    return axis - 1


def cell_centroids(mesh):
    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim, 0)
    c_to_v = mesh.topology.connectivity(tdim, 0)
    x = mesh.geometry.x
    cell_map = mesh.topology.index_map(tdim)
    n_cells = cell_map.size_local + cell_map.num_ghosts
    centroids = np.zeros((n_cells, 3), dtype=float)
    for cell in range(n_cells):
        vertices = c_to_v.links(cell)
        centroids[cell, :] = x[vertices, :3].mean(axis=0)
    return centroids


def cell_global_indices(mesh):
    tdim = mesh.topology.dim
    cell_map = mesh.topology.index_map(tdim)
    n_cells = cell_map.size_local + cell_map.num_ghosts
    local_cells = np.arange(n_cells, dtype=np.int32)
    return np.asarray(cell_map.local_to_global(local_cells), dtype=np.int64)


def global_bbox_center(points: np.ndarray):
    local_min = points.min(axis=0)
    local_max = points.max(axis=0)
    global_min = np.empty(3, dtype=float)
    global_max = np.empty(3, dtype=float)
    comm.Allreduce(local_min, global_min, op=MPI.MIN)
    comm.Allreduce(local_max, global_max, op=MPI.MAX)
    return 0.5 * (global_min + global_max)


def rve_centroid_inplane(fiber_axis: int):
    center = global_bbox_center(domain.geometry.x[:, :3])
    axis = _axis_index(fiber_axis)
    if axis == 0:
        return center[2], center[1]  # Fortran: a0=z, b0=y
    if axis == 1:
        return center[0], center[2]  # Fortran: a0=x, b0=z
    return center[0], center[1]  # Fortran: a0=x, b0=y


def owned_partition_payload(global_cells, centroids, matrix_cells, fiber_cells, coord_axis):
    n_owned = domain.topology.index_map(domain.topology.dim).size_local
    matrix_owned = matrix_cells[matrix_cells < n_owned]
    fiber_owned = fiber_cells[fiber_cells < n_owned]

    fiber_payload = np.column_stack(
        (global_cells[fiber_owned], np.zeros(len(fiber_owned), dtype=float))
    )
    matrix_payload = np.column_stack(
        (global_cells[matrix_owned], centroids[matrix_owned, coord_axis])
    )
    return fiber_payload, matrix_payload


def gather_payload(local_payload):
    gathered = comm.allgather(local_payload)
    nonempty = [a for a in gathered if a.size]
    if not nonempty:
        return np.empty((0, 2), dtype=float)
    return np.vstack(nonempty)


def part_id_from_global_map(global_cells, global_part_map):
    part_id = np.zeros(len(global_cells), dtype=np.int32)
    for local_cell, global_cell in enumerate(global_cells):
        part_id[local_cell] = global_part_map.get(int(global_cell), 0)
    return part_id


def assign_strip_partitions(centroids, global_cells, matrix_cells, fiber_cells, m_parts, strip_axis):
    if m_parts < 1:
        raise ValueError("M_PARTS must be >= 1")

    fiber_payload, matrix_payload = owned_partition_payload(
        global_cells, centroids, matrix_cells, fiber_cells, _axis_index(strip_axis)
    )
    all_fiber = gather_payload(fiber_payload)
    all_matrix = gather_payload(matrix_payload)
    if m_parts > len(all_matrix):
        raise ValueError("M_PARTS cannot exceed the global number of matrix cells")

    order = np.lexsort((all_matrix[:, 0], all_matrix[:, 1]))
    sorted_global_matrix = all_matrix[order, 0].astype(np.int64)
    global_part_map = {int(g): 1 for g in all_fiber[:, 0].astype(np.int64)}
    for strip, global_ids in enumerate(np.array_split(sorted_global_matrix, m_parts), start=2):
        for global_id in global_ids:
            global_part_map[int(global_id)] = strip
    return part_id_from_global_map(global_cells, global_part_map), m_parts + 1


def assign_angular_partitions(centroids, global_cells, matrix_cells, fiber_cells, m_parts, fiber_axis):
    if m_parts < 1:
        raise ValueError("M_PARTS must be >= 1")

    a0, b0 = rve_centroid_inplane(fiber_axis)
    axis = _axis_index(fiber_axis)
    if axis == 0:
        ca = centroids[matrix_cells, 1] - b0
        cb = centroids[matrix_cells, 2] - a0
    elif axis == 1:
        ca = centroids[matrix_cells, 0] - a0
        cb = centroids[matrix_cells, 2] - b0
    else:
        ca = centroids[matrix_cells, 0] - a0
        cb = centroids[matrix_cells, 1] - b0
    theta = np.arctan2(cb, ca)
    theta = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)
    sector = np.floor(theta / (2.0 * np.pi / m_parts)).astype(np.int32) + 2
    sector = np.clip(sector, 2, m_parts + 1)

    part_id = np.zeros(len(global_cells), dtype=np.int32)
    part_id[fiber_cells] = 1
    part_id[matrix_cells] = sector
    return part_id, m_parts + 1

def square_radius(theta, half_x=0.5, half_y=0.5):
    c = np.abs(np.cos(theta))
    s = np.abs(np.sin(theta))
    eps = 1e-14
    return np.minimum(
        half_x / np.maximum(c, eps),
        half_y / np.maximum(s, eps),
    )


def assign_radial_angular_partitions(
    centroids, global_cells, matrix_cells, fiber_cells,
    n_angular, n_radial, fiber_axis
):
    a0, b0 = rve_centroid_inplane(fiber_axis)
    axis = _axis_index(fiber_axis)

    if axis == 0:
        ca = centroids[matrix_cells, 1] - b0
        cb = centroids[matrix_cells, 2] - a0
        half_a, half_b = 0.5, 0.5
    elif axis == 1:
        ca = centroids[matrix_cells, 0] - a0
        cb = centroids[matrix_cells, 2] - b0
        half_a, half_b = 0.5, 0.5
    else:
        ca = centroids[matrix_cells, 0] - a0
        cb = centroids[matrix_cells, 1] - b0
        half_a, half_b = 0.5, 0.5

    theta = np.arctan2(cb, ca)
    theta = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)

    r = np.sqrt(ca**2 + cb**2)

    r_inner = r.min()
    r_outer = square_radius(theta, half_a, half_b)

    rho = (r - r_inner) / (r_outer - r_inner + 1e-14)
    rho = np.clip(rho, 0.0, 1.0)

    angular_id = np.floor(theta / (2.0 * np.pi / n_angular)).astype(np.int32)
    angular_id = np.clip(angular_id, 0, n_angular - 1)

    radial_id = np.floor(rho * n_radial).astype(np.int32)
    radial_id = np.clip(radial_id, 0, n_radial - 1)

    # matrix_part = radial_id * n_angular + angular_id + 2
    matrix_part = angular_id * n_radial + radial_id + 2
    part_id = np.zeros(len(global_cells), dtype=np.int32)
    part_id[fiber_cells] = 1
    part_id[matrix_cells] = matrix_part

    n_parts = 1 + n_angular * n_radial
    return part_id, n_parts


def part_indicator(part_id, part_number, name):
    chi = fem.Function(W)
    chi.name = name
    chi.x.array[:] = 0.0
    chi.x.array[np.flatnonzero(part_id == part_number)] = 1.0
    return chi


def collect_partition_rows(global_cells, part_id):
    if rank != 0:
        return None
    n_owned = domain.topology.index_map(domain.topology.dim).size_local
    rows = np.column_stack((global_cells[:n_owned], part_id[:n_owned]))
    return rows[np.argsort(rows[:, 0])].astype(np.int64)


def write_partition_file(filename, partition_rows):
    if rank != 0 or not filename:
        return
    with open(filename, "w") as f:
        f.write("# global_cell_id   part_id\n")
        for global_cell, pid in partition_rows:
            f.write(f"{global_cell + 1:14d}  {pid:5d}\n")


def write_area_file(filename, volumes):
    if rank != 0 or not filename:
        return
    vol_frac = np.asarray(volumes, dtype=float)
    total = vol_frac.sum()
    if total <= 0.0:
        raise RuntimeError("Cannot write area.dat: non-positive total partition volume")
    vol_frac = vol_frac / total
    with open(filename, "w") as f:
        f.write(" ".join(f"{v:.16g}" for v in vol_frac))
        f.write("\n")


centroids = cell_centroids(domain)
global_cells = cell_global_indices(domain)


if PART_TYPE.lower() == "radial_angular":
    part_id, n_parts = assign_radial_angular_partitions(
        centroids, global_cells, mat_cells, fib_cells,
        N_ANGULAR, N_RADIAL, FIBER_AXIS
    )
elif PART_TYPE.lower() == "angular":
    part_id, n_parts = assign_angular_partitions(
        centroids, global_cells, mat_cells, fib_cells, M_PARTS, FIBER_AXIS
    )
else:
    part_id, n_parts = assign_strip_partitions(
        centroids, global_cells, mat_cells, fib_cells, M_PARTS, STRIP_AXIS
    )

if np.any(part_id == 0):
    raise RuntimeError("Some local cells were not assigned to an EBH partition")

partition_names = ["fiber"] + [f"matrix_{PART_TYPE.lower()}_{i}" for i in range(1, M_PARTS + 1)]
partitions = [
    (partition_names[p - 1], part_indicator(part_id, p, f"chi_{partition_names[p - 1]}"))
    for p in range(1, n_parts + 1)
]
partition_rows = collect_partition_rows(global_cells, part_id)
write_partition_file(PART_FILE, partition_rows)

part_id_fun = fem.Function(W)
part_id_fun.name = "partition_id"
part_id_fun.x.array[:] = part_id.astype(float)
if rank == 0:
    with XDMFFile(comm, "vtk_outputs/partitions.xdmf", "w") as xdmf_part:
        xdmf_part.write_mesh(domain)
        xdmf_part.write_function(part_id_fun)

# Volumes of partitions
volumes = []
dx = ufl.dx(domain=domain)
if rank == 0:
    print(
        f"\nPartitioning: PART_TYPE={PART_TYPE}, M_PARTS={M_PARTS}, "
        f"n_parts={n_parts}"
    )
n_owned_cells = domain.topology.index_map(domain.topology.dim).size_local
for p, (name, chi) in enumerate(partitions, start=1):
    local_vol = fem.assemble_scalar(fem.form(chi * 1.0 * dx))
    vol = comm.allreduce(local_vol, op=MPI.SUM)
    volumes.append(vol)
    n_local = int(np.count_nonzero(part_id[:n_owned_cells] == p))
    n_global = n_local
    # if rank == 0:
    #     print(f"partition {p:3d} = {name:>18s}, cells = {n_global:7d}, volume = {vol:.8e}")


write_area_file(AREA_FILE, volumes)
if rank == 0:
    print(f"  written: {AREA_FILE}")


# -----------------------------------------------------------------------------
# Linear elasticity bilinear form
# -----------------------------------------------------------------------------
du = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a_form = ufl.inner(C_eps(eps(du)), eps(v)) * dx

petsc_options = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}
# Each rank solves serial local systems in task-parallel mode, so no parallel LU
# package is required.


class ReusableMPCLinearSolver:
    """Assemble/factorize the MPC stiffness matrix once, reuse it for many RHS."""

    def __init__(self, a_ufl, mpc, bcs, petsc_options):
        self.a = fem.form(a_ufl)
        self.mpc = mpc
        self.bcs = bcs
        self.A = assemble_matrix(self.a, self.mpc, bcs=self.bcs)
        self.A.assemble()

        self.b = self.A.createVecRight()
        self.x = self.A.createVecRight()
        self.u = fem.Function(self.mpc.function_space)

        self.ksp = PETSc.KSP().create(self.A.comm)
        self.ksp.setOperators(self.A)
        prefix = f"ebh_reuse_{id(self)}_"
        self.ksp.setOptionsPrefix(prefix)
        opts = PETSc.Options()
        opts.prefixPush(prefix)
        for key, value in petsc_options.items():
            opts[key] = value
        self.ksp.setFromOptions()
        for key in petsc_options:
            del opts[key]
        opts.prefixPop()

        self.ksp.setUp()
        pc = self.ksp.getPC()
        try:
            factor_solver = pc.getFactorSolverType()
        except PETSc.Error:
            factor_solver = "not-a-factor-pc"
        print(
            f"[rank {rank:02d}] PETSc solver: "
            f"ksp={self.ksp.getType()}, pc={pc.getType()}, factor={factor_solver}",
            flush=True,
        )

    def solve(self, L):
        self.b.set(0.0)
        assemble_vector(L, self.mpc, b=self.b)
        apply_lifting(self.b, [self.a], bcs=[self.bcs], constraint=self.mpc)
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(self.b, self.bcs)
        self.b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        self.x.set(0.0)
        self.ksp.solve(self.b, self.x)
        self.x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        with self.x.localForm() as x_local:
            self.u.x.array[: x_local.array_r.size] = x_local.array_r
        self.u.x.scatter_forward()
        self.mpc.homogenize(self.u)
        self.mpc.backsubstitution(self.u)

        u_out = fem.Function(self.mpc.function_space)
        u_out.x.array[:] = self.u.x.array
        u_out.x.scatter_forward()
        return u_out


solver = ReusableMPCLinearSolver(a_form, mpc, [bc_pin], petsc_options)
if rank == 0:
    print("\nAssembled and factorized the MPC stiffness matrix once per rank.", flush=True)


# -----------------------------------------------------------------------------
# Precompiled RHS and averaging forms
# -----------------------------------------------------------------------------
E_units = [unit_E(kl) for kl in range(NV)]
E_rhs_forms = [
    fem.form(-ufl.inner(C_eps(ufl.as_tensor(E_units[kl])), eps(v)) * dx)
    for kl in range(NV)
]

chi_rhs = fem.Function(W)
chi_rhs.name = "chi_rhs_source_partition"
P_rhs_forms = [
    fem.form(
        ufl.inner(C_eps(masked_tensor(E_units[kl], chi_rhs)), eps(v)) * dx
    )
    for kl in range(NV)
]

avg_u = fem.Function(V)
avg_strain_forms = [
    [
        [fem.form(chi_B * eps(avg_u)[i, j] * dx) for j in range(3)]
        for i in range(3)
    ]
    for _, chi_B in partitions
]


def average_strain_precomputed(u_sol, B, add_macro=None):
    avg_u.x.array[:] = u_sol.x.array
    avg_u.x.scatter_forward()

    M = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            local_val = fem.assemble_scalar(avg_strain_forms[B][i][j])
            M[i, j] = comm.allreduce(local_val, op=MPI.SUM) / volumes[B]
    if add_macro is not None:
        M += add_macro
    return to_voigt(M)


if rank == 0:
    print("Precompiled RHS forms and partition-average strain forms.", flush=True)

world_comm.Barrier()
t_parallel_solve_start = MPI.Wtime()


# -----------------------------------------------------------------------------
# Compute E_inf and P_inf
# -----------------------------------------------------------------------------
E_inf = np.zeros((NV, NV, n_parts), dtype=float)
P_inf = np.zeros((NV, NV, n_parts, n_parts), dtype=float)

H_modes = []
h_fields = [[None for _ in range(NV)] for _ in range(n_parts)]

xdmf_H = None
xdmf_P = None
if rank == 0:
    xdmf_H = XDMFFile(comm, "vtk_outputs/E_H_modes.xdmf", "w")
    xdmf_H.write_mesh(domain)
    xdmf_P = XDMFFile(comm, "vtk_outputs/P_h_modes_rank0.xdmf", "w")
    xdmf_P.write_mesh(domain)

local_e_task_count = sum(1 for kl in range(NV) if kl % size == rank)
local_e_done_count = 0
all_e_task_counts = world_comm.gather(local_e_task_count, root=0)
if rank == 0:
    print("\nSolving macro-strain problems for E_inf...", flush=True)
    print("  E tasks per rank:", all_e_task_counts, flush=True)
world_comm.Barrier()

for kl in range(NV):
    if kl % size != rank:
        continue

    Eunit = E_units[kl]
    # Fortran equivalent: K u = -F_macro(:, kl)
    uH = solver.solve(E_rhs_forms[kl])

    uH.name = f"H_macro_mode_{VOIGT_LABELS[kl]}"
    # if rank == 0:
    #     xdmf_H.write_function(uH)

    H_modes.append(uH)
    for B, (_, chi_B) in enumerate(partitions):
        # E_inf includes the macro identity term: <eps(H) + Eunit>_B
        E_inf[:, kl, B] = average_strain_precomputed(
            uH, B, add_macro=Eunit
        )

    local_e_done_count += 1
    print(
        f"[rank {rank:02d}] finished E task {local_e_done_count}/{local_e_task_count}: "
        f"mode={VOIGT_LABELS[kl]}",
        flush=True,
    )

E_global = np.zeros_like(E_inf) if rank == 0 else None
world_comm.Reduce(E_inf, E_global, op=MPI.SUM, root=0)
if rank == 0:
    E_inf = E_global

total_p_tasks = n_parts * NV
local_task_count = sum(1 for task_id in range(total_p_tasks) if task_id % size == rank)
local_done_count = 0
all_task_counts = world_comm.gather(local_task_count, root=0)
if rank == 0:
    print("\nSolving partition eigenstrain problems for P_inf...", flush=True)
    print("  task-parallel ranks:", size, flush=True)
    print("  P tasks per rank:", all_task_counts, flush=True)
world_comm.Barrier()

for A, (part_name_A, chi_A) in enumerate(partitions):
    for kl in range(NV):
        task_id = A * NV + kl
        if task_id % size != rank:
            continue

        print(
            f"[rank {rank:02d}] starting P task {local_done_count + 1}/{local_task_count}: "
            f"A={A + 1}/{n_parts} ({part_name_A}), mode={VOIGT_LABELS[kl]}",
            flush=True,
        )

        chi_rhs.x.array[:] = chi_A.x.array
        chi_rhs.x.scatter_forward()

        # Fortran equivalent: K u = F_eig(A, kl)
        h_sol = solver.solve(P_rhs_forms[kl])
        h_sol.name = f"h_A{A + 1}_{part_name_A}_mode_{VOIGT_LABELS[kl]}"
        # if rank == 0:
        #     xdmf_P.write_function(h_sol)
        h_fields[A][kl] = h_sol

        for B, (_, chi_B) in enumerate(partitions):
            # P_inf is only the fluctuation strain average <eps(h)>_B
            P_inf[:, kl, B, A] = average_strain_precomputed(
                h_sol, B, add_macro=None
            )

        local_done_count += 1
        print(
            f"[rank {rank:02d}] finished P task {local_done_count}/{local_task_count}: "
            f"A={A + 1}/{n_parts} ({part_name_A}), mode={VOIGT_LABELS[kl]}",
            flush=True,
        )

P_global = np.zeros_like(P_inf) if rank == 0 else None
world_comm.Reduce(P_inf, P_global, op=MPI.SUM, root=0)
if rank == 0:
    P_inf = P_global

world_comm.Barrier()
elapsed_parallel_solve = world_comm.reduce(
    MPI.Wtime() - t_parallel_solve_start, op=MPI.MAX, root=0
)

if rank == 0:
    xdmf_H.close()
    xdmf_P.close()


# -----------------------------------------------------------------------------
# Write outputs with the same layout as the Fortran code
# -----------------------------------------------------------------------------
def write_E_inf(filename, E_inf):
    if rank != 0:
        return
    with open(filename, "w") as f:
        f.write("# E_inf(ij, kl, B) -- rows: (B-1)*6+ij, cols: kl\n")
        f.write("# partition order: " + ", ".join(name for name, _ in partitions) + "\n")
        f.write("# voigt order: xx yy zz xy yz xz\n")
        for B in range(n_parts):
            for ij in range(NV):
                f.write(" ".join(f"{E_inf[ij, kl, B]: .14E}" for kl in range(NV)))
                f.write("\n")


def write_P_inf(filename, P_inf):
    if rank != 0:
        return
    with open(filename, "w") as f:
        f.write("# P_inf(ij, kl, B, A) -- rows: (B-1)*6+ij, cols: (A-1)*6+kl\n")
        f.write("# partition order: " + ", ".join(name for name, _ in partitions) + "\n")
        f.write("# voigt order: xx yy zz xy yz xz\n")
        for B in range(n_parts):
            for ij in range(NV):
                row = []
                for A in range(n_parts):
                    for kl in range(NV):
                        row.append(f"{P_inf[ij, kl, B, A]: .14E}")
                f.write(" ".join(row))
                f.write("\n")


write_E_inf("E_inf.dat", E_inf)
write_P_inf("P_inf.dat", P_inf)

if rank == 0:
    np.savez(
        "roh_cache/ebh_offline_E_P.npz",
        E_inf=E_inf,
        P_inf=P_inf,
        partition_names=np.array(partition_names),
        part_id=partition_rows[:, 1],
        global_cell_id=partition_rows[:, 0],
        part_type=np.array(PART_TYPE),
        matrix_partition_count=np.array(M_PARTS),
        partition_volumes=np.array(volumes),
        partition_volume_fractions=np.array(volumes) / np.sum(volumes),
        voigt_order=VOIGT_LABELS,
    )

    print("\nEBH offline done.")
    print("  written: E_inf.dat")
    print("  written: P_inf.dat")
    print(f"  written: {PART_FILE}")
    print(f"  written: {AREA_FILE}")
    print("  written: vtk_outputs/partitions.xdmf")
    print("  written: roh_cache/ebh_offline_E_P.npz")
    print("  E_inf shape =", E_inf.shape)
    print("  P_inf shape =", P_inf.shape)

world_comm.Barrier()
elapsed_total = world_comm.reduce(MPI.Wtime() - t_total_start, op=MPI.MAX, root=0)
if rank == 0:
    print(f"  parallel solve wall time (no assembly/precompute/write) = {elapsed_parallel_solve:.3f} s")
    print(f"  total wall time = {elapsed_total:.3f} s")
