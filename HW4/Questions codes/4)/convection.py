from firedrake import *
import firedrake_ts
import os

os.makedirs("result", exist_ok=True)

# -------------------
# Parameters
# -------------------
N = 128
Ra = Constant(1.0e4)
A = Constant(0.1)

t_init = 0.0
t_max = 100000.0

mesh = UnitSquareMesh(N, N, quadrilateral=True)

V = FunctionSpace(mesh, "CG", 1)
ME = MixedFunctionSpace([V, V, V], name=["temperature", "vorticity", "streamfunction"])
Vv = VectorFunctionSpace(mesh, "CG", 1)

# Unknowns and time derivatives
u = Function(ME, name="solution")
u_dot = Function(ME, name="solution_dot")

T, omega, psi = split(u)
T_dot, omega_dot, psi_dot = split(u_dot)

q, eta, phi = TestFunctions(ME)

x, z = SpatialCoordinate(mesh)

# -------------------
# Initial condition
# -------------------
T0 = 1.0 - z + A*cos(pi*x)#*sin(pi*x)

u.subfunctions[0].interpolate(T0)
u.subfunctions[1].interpolate(0.0)
u.subfunctions[2].interpolate(0.0)

u.subfunctions[0].rename("Temperature")
u.subfunctions[1].rename("Vorticity")
u.subfunctions[2].rename("Streamfunction")

# Velocity: v = curl(psi k)
vel = as_vector((psi.dx(1), -psi.dx(0)))

# -------------------
# Weak form
# -------------------

# Temperature:
# T_t + v dot grad(T) - (1/Ra) Delta T = 0
FT = (q*T_dot*dx+ q*dot(vel, grad(T))*dx + (1.0/Ra)*inner(grad(T), grad(q))*dx)
Fomega = (inner(grad(eta),grad(omega))*dx-eta*T.dx(0)*dx)
Fpsi = (inner(grad(phi),grad(psi))*dx-phi*omega*dx)
F = FT + Fomega + Fpsi

# -------------------
# Boundary conditions
# -------------------
# Firedrake unit square boundary IDs:
# 1: x=0, 2: x=1, 3: y=0, 4: y=1

bcs = [
    DirichletBC(ME.sub(0), 1.0, 3),   # T(x,0,t)=1
    DirichletBC(ME.sub(0), 0.0, 4),   # T(x,1,t)=0
    DirichletBC(ME.sub(1), 0.0, "on_boundary"),  # omega=0
    DirichletBC(ME.sub(2), 0.0, "on_boundary"),  # psi=0
]


# -------------------
# Output
# -------------------
outfile = VTKFile(f"result/convection_Ra_{float(Ra):.0e}_N_{N}.pvd")

velocity_out = Function(Vv, name="Velocity")

Nu_file = open(f"result/nusselt_Ra_{float(Ra):.0e}_N_{N}.txt", "w")

def compute_nusselt(Tfun):
    h = 1.0  # altura del dominio
    numerator   = assemble(Tfun.dx(1) * ds(4))  # ∂T/∂z en z=1
    denominator = assemble(Tfun * ds(3))         # T en z=0 (debería ser ≈1)
    return -h*numerator/denominator

def monitor(ts, step, t, x_petsc):
    Tfun = u.subfunctions[0]
    psifun = u.subfunctions[2]

    # velocity_out.interpolate(curl(psifun))
    velocity_out.interpolate(as_vector((psifun.dx(1), -psifun.dx(0))))
    if step % 10 == 0:
        Nu = compute_nusselt(Tfun)
        print(f"step={step}, t={t:.6e}, Nu={Nu:.8f}")
        Nu_file.write(f"{step} {t:.16e} {Nu:.16e}\n")
        Nu_file.flush()

        outfile.write(
            u.subfunctions[0],
            u.subfunctions[1],
            u.subfunctions[2],
            velocity_out,
            time=t,
        )

# -------------------
# Solver parameters
# -------------------
params = {
    "ts_type": "bdf",
    "ts_bdf_order": 2,
    "ts_dt": 1,
    # "ts_monitor": None,
    "ts_rtol": 1e-6,
    "ts_atol": 1e-10,
    "ts_max_time": t_max,
    "ts_adapt_dt_min": 1e-3,
    "ts_adapt_dt_max": 100.0,
    "ts_exact_final_time": "matchstep",
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}

problem = firedrake_ts.DAEProblem(
    F,
    u,
    u_dot,
    (t_init, t_max),
    bcs=bcs,
)

solver = firedrake_ts.DAESolver(
    problem,
    solver_parameters=params,
    monitor_callback=monitor,
)

solver.solve()

Nu_final = compute_nusselt(u.subfunctions[0])
print(f"\nFinal Nu = {Nu_final:.10f}")

Nu_file.close()