# Solve the biharmonic equation as a coupled system of two Poisson equations.
# This code reproduces the manufactured solution for the biharmonic equation from the Finite Difference code in C,
#  but using the Firedrake finite element library.
import firedrake_ts
from firedrake import *
import numpy as np


N = 2 
levels = 8
Nfine = N*2**levels
base_mesh = UnitSquareMesh(N, N, quadrilateral=True)
Hierarchy = MeshHierarchy(base_mesh, levels)
mesh = Hierarchy[-1]
V = FunctionSpace(mesh, "Lagrange", 1)
ME = MixedFunctionSpace([V, V], name=["vorticity", "streamfunction"])
Vv = VectorFunctionSpace(mesh, "Lagrange", 1)


# Define test functions for porosity and pressure
omega_t, psi_t = TestFunctions(ME)

# Define functions
u = Function(ME)  # current solution
u.subfunctions[0].rename("vorticity")
u.subfunctions[1].rename("streamfunction")

# Split mixed functions
omega, psi = split(u)

# set RHS function and True solutions
f = Function(V, name="rhs")
# set initial conditions
x, z = SpatialCoordinate(mesh)
# Define the manufactured solution and corresponding RHS for the biharmonic equation


A = 0.1
T = Function(V, name="Temperature")
T.interpolate((1.0 - z) + A * cos(pi * x))
f.interpolate(T.dx(0))

# Weak statement of the equations

Fomega = inner( grad(omega_t), (grad(omega)) ) * dx - omega_t * f * dx
Fpsi = inner( grad(psi_t), (grad(psi)) ) * dx -  psi_t * omega * dx
F = Fomega + Fpsi


# set initial conditions
u.subfunctions[0].interpolate(0.)
u.subfunctions[1].interpolate(0.)

# boundary conditions
bcs = [DirichletBC(ME.sub(0), 0., "on_boundary"),
       DirichletBC(ME.sub(1), 0., "on_boundary")]

# pc = "monolithic" 
pc = "fieldsplit_direct"
# pc = "fieldsplit_mg"
params_general = {
    "snes_type": "ksponly",
    "ksp_rtol": 1e-6,
    "ksp_atol": 1e-10,
    "snes_monitor": None,
    "ksp_monitor": None
}  
if pc == "monolithic": 
    params = {  
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps"
    } 
elif pc == "fieldsplit_mg":
    params = {
        "ksp_type": "fgmres",
        "pc_type": "fieldsplit",
        "pc_fieldsplit_type": "multiplicative",  # Options: additive, multiplicative, schur
        # Nested dictionaries for solver settings on each block
        "fieldsplit_0": {
            "ksp_type": "preonly",
            "pc_type": "mg",
            "pc_factor_mat_solver_type": "mumps"
        },
        "fieldsplit_1": {
            "ksp_type": "preonly",
            "pc_type": "mg",
            "pc_factor_mat_solver_type": "mumps"
        }
    } 
elif pc == "fieldsplit_direct":
    params = {
        "ksp_type": "fgmres",
        "pc_type": "fieldsplit",
        "pc_fieldsplit_type": "multiplicative",
        "fieldsplit_0": {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps"
        },
        "fieldsplit_1": {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps"
        }
    }

params.update(params_general)
#print(params)

problem = NonlinearVariationalProblem(F, u, bcs=bcs)
solver = NonlinearVariationalSolver(problem, solver_parameters=params)

#print(dir(problem))
solver.solve()

v = Function(Vv, name="Velocity")
v.interpolate(curl(psi))
outfile = VTKFile("result/biharm.pvd")
outfile.write(T, f, u.subfunctions[0], u.subfunctions[1], v)

#calculate error norms
print(f"\nworking on {Nfine}x{Nfine} mesh with {pc} preconditioner and {params['ksp_type']} solver")


