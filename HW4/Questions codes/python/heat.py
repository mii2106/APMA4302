from firedrake import *
import firedrake_ts

N=4 
levels=4
coarsemesh = UnitSquareMesh(N, N, quadrilateral=True)
hierarchy = MeshHierarchy(coarsemesh, levels)
mesh = hierarchy[-1]




family = "CG"
degree = 1      
V = FunctionSpace(mesh, family, degree)

u = Function(V, name="Temperature")
u_dot  = Function(V, name="u_dot")
v = TestFunction(V)

# set rhs source function to constant
x, y = SpatialCoordinate(mesh)
f = Function(V)
f.interpolate(1.0)

F = inner(u_dot, v) * dx + inner(grad(u), grad(v)) * dx - f * v * dx

# set Dirichlet boundary conditions on left and bottom boundaries
bc_val = Constant(0.0)
bc = DirichletBC(V, bc_val, (1, 3))

t_init = 0.0
t_max = 4.0
params = {'ts_type': 'bdf',
          'ts_bdf_order': 2,
          'ts_dt': 0.01,
          'ts_monitor': None,
          'ts_rtol': 1e-6,
          'ts_atol': 1e-10,
          'ksp_type': 'preonly',
          'pc_type': 'lu',
          'pc_factor_mat_solver_type': 'mumps',
          'ts_max_time': t_max,
          'ts_adapt_dt_min': 1.e-6,
          'ts_exact_final_time': "matchstep"}

outfile = VTKFile("result/heat.pvd")
outfile.write(u, time=0.0)

def monitor(ts, step, t, x):
    outfile.write(u, time=t)    

problem = firedrake_ts.DAEProblem(F, u, u_dot, (t_init, t_max), bcs=bc)
solver = firedrake_ts.DAESolver(problem, solver_parameters=params, monitor_callback=monitor)

solver.solve()
