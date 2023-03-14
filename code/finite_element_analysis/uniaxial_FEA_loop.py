##################
#
# Uniaxial stretch FEA code
# by M.S.Dimitriyev
# michael.dimitriyev@gmail.com
#
# Uses FEniCS 2019.1.0 as
# packaged in Anaconda
#
##################
from dolfin import *
import numpy as np
import sys
import subprocess


# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True,
               "eliminate_zeros": True, 
               "precompute_basis_const": True,
               "precompute_ip_const": True}

# 40.25 mm tall, 126.66 mm wide
L = 3.15
H = 1.0
Nx = 50
Ny = 50

# Max displacement and displacement step for sweep
DEF_MAX = 1.2
dy = 0.02


for i in range(0,int(DEF_MAX/dy)+1):
	DEF_RATIO = round(i*dy,2)
	DISP = DEF_RATIO*H
	
	out = "out_" + str(DEF_RATIO)
	subprocess.run(["mkdir",out])

	mesh = RectangleMesh(Point(0., 0.), Point(L, H), Nx, Ny, "crossed")
	File("my_mesh.xml") << mesh
	
	# Create mesh and define function space
	V = VectorFunctionSpace(mesh, "Lagrange", 1)
	
	# Create boundary subdomains
	bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)  
	top = CompiledSubDomain("near(x[1], side) && on_boundary", side = H)
	
	# Define Dirichlet boundary (y = 0, y = 1)
	b = Expression(("0.0", "0.0"),degree=1)
	t = Expression(("0.0", "disp"),disp=DISP,degree=1)
	
	# Assign boundary conditions
	bcb = DirichletBC(V, b, bottom)
	bct = DirichletBC(V, t, top)
	bcs = [bcb, bct]
	
	
	# Define functions
	du = TrialFunction(V)            # Incremental displacement
	v  = TestFunction(V)             # Test function
	u  = Function(V)                 # Displacement from previous iteration
	B  = Constant((0.0, 0.0))  # Body force per unit volume
	T  = Constant((0.0, 0.0))  # Traction force on the boundary
	
	# Define strain
	d = u.geometric_dimension()
	strain = 0.5*(grad(u) + grad(u).T)
	
	# use Voigt notation to enable matrix multiplication for anisotropy
	e_vgt = as_vector([strain[0,0],strain[1,1],2*strain[0,1]])
	

	# anisotropic elasticity tensor
	e1 = Constant( (1.0, 0.0, 0.0) )
	e2 = Constant( (0.0, 1.0, 0.0) )
	e3 = Constant( (0.0, 0.0, 1.0) )
	c_hook1  = Constant( ( (0.240872, 0.0362245, 0.0), (0.0362245, 0.0599528, 0.0), (0.0, 0.0, 0.01) ) )
	alpha_xx = Constant( (1.16995, 0.0, 0.0) )
	alpha_yy = Constant( (0.0, 0.802057, 0.0) )
	beta_xx  = Constant( (0.0223765, 0.0, 0.0) )
	beta_yy  = Constant( (0.0, 0.0219965, 0.0) )	


	# Stored strain energy density
	psi = .5*dot(e_vgt,c_hook1*e_vgt) 
	psi += dot(e1,beta_xx)*( pow(dot(e1,alpha_xx),2)*pow(dot(e1,e_vgt),3) + pow(dot(e1,alpha_xx),3)*pow(dot(e1,e_vgt),4) )
	psi += dot(e2,beta_yy)*( pow(dot(e2,alpha_yy),2)*pow(dot(e2,e_vgt),3) + pow(dot(e2,alpha_yy),3)*pow(dot(e2,e_vgt),4) )

	# Total potential energy
	Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
	
	# Compute first variation of Pi (directional derivative about u in the direction of v)
	F = derivative(Pi, u, v)
	
	# Compute Jacobian of F
	J = derivative(F, u, du)
	
	# Solve variational problem
	solve(F == 0, u, bcs, J=J,
	      form_compiler_parameters=ffc_options,
	      solver_parameters={
	      "newton_solver":{"relative_tolerance":1E-12},
	      "newton_solver":{"maximum_iterations":100}})
	
	
	# Write to file
	X = V.tabulate_dof_coordinates()
	X.resize((V.dim(), 2))
	
	x = X[:,0]
	y = X[:,1]
	
	disp_out = open(out+"/disp.out","w")
	disp_field  = u.vector().get_local()
	for it in range(0,int(disp_field.size/2)):
		disp_out.write(str(x[2*it]) + ", " + str(y[2*it]) + ", " + str(disp_field[2*it]) + ", " + str(disp_field[2*it+1]) + "\n")
	disp_out.close()
	
	scalar_ed = FunctionSpace(mesh, 'DG', 0)
	
	ed_field = project(psi, scalar_ed).vector().get_local()
	
	ed_out = open(out+"/ed.out","w")
	for it in range(0,int(ed_field.size)):
		ed_out.write(str(ed_field[it]) + "\n")
	ed_out.close()
	

	energy_out = open(out+"/energy.out","w")
	energy_out.write(str(DEF_RATIO) + "\t" + str(assemble(Pi)))
	energy_out.close()
