import coolfluid as cf
import math

### Create new model specialized for SD
model   = cf.root.create_component('cylinder_2d','cf3.dcm.Model');

### Load the mesh
# mesh = model.domain.load_mesh(file = cf.URI('../../../resources/cylinder-quad-p2-16x4.msh'), name = 'cylinder2d');
mesh = model.domain.load_mesh(file = cf.URI('../../../resources/cylinder-quad-p2-32x8.msh'), name = 'cylinder2d');
# mesh = model.domain.load_mesh(file = cf.URI('../../../resources/cylinder-quad-p2-64x16.msh'), name = 'cylinder2d');
# mesh = model.domain.load_mesh(file = cf.URI('../../../resources/cylinder-quad-p2-128x32.msh'), name = 'cylinder2d');
model.build_faces();

### Compute some physics variables
gamma = 1.4
R = 287.05
M_inf = 0.38
p_inf = 1.0
rho_inf = 1.0
c_inf = math.sqrt(gamma*p_inf/rho_inf)
u_inf = M_inf*c_inf
rhoE_inf = p_inf/(gamma-1) + 0.5 * rho_inf * u_inf**2
#p_inf = (g-1) * ( rhoE - 0.5 * rho * ( u**2 ) )

### Add the Partial Differential Equations to solve
euler = model.add_pde(name='euler',type='cf3.dcm.equations.euler.Euler2D',shape_function='cf3.dcm.core.LegendreGaussEndP1')
euler.gamma = gamma
euler.R = R

euler.add_term(name='rhs',type='cf3.sdm.br2.euler_RightHandSide2D')

### Add BC
bc_wall = euler.add_bc( name='wall',
                        type='cf3.dcm.equations.euler.BCMirror2D',
                        regions=[ mesh.topology.cylinder ] )
bc_farfield = euler.add_bc( name = 'farfield',
                            type = 'cf3.dcm.equations.euler.BCFarfield2D',
                            regions= [ mesh.topology.boundary ] )
bc_farfield.rho = rho_inf
bc_farfield.U   = [u_inf,0]
bc_farfield.p   = p_inf


### Initialize the solution
model.tools.init_field.init_field(
  field=euler.solution,
  functions=[
    str(rho_inf),
    'if( (x<-5) | (x>5) , '+str(rho_inf*u_inf)+' , 0.5*'+str(rho_inf*u_inf)+' )',
    '0.',
    str(rhoE_inf)] )

### Create the Solver for the Partial Differential Equations
solver = model.add_solver(name='rk_solver',pde=euler)
solver.children.time_step_computer.cfl = 0.5
solver.options.order = 3


solution  = euler.fields.solution
post_proc = euler.fields.create_field(name='post_proc',
                                      variables='U[vec],p[1],T[1],M[1],Pt[1],Tt[1],Cp[1],S[1]')

### Time Stepping
if True:
  # solver.solve_iterations(500)
  solver.solve_time_step(1)
  
  for index in range(len(solution)) :
     # compute variables per solution point
     rho=solution[index][0];
     u=solution[index][1]/solution[index][0];
     v=solution[index][2]/solution[index][0];
     rhoE=solution[index][3];
     p=(gamma-1)*(rhoE - 0.5*rho*(u**2+v**2));
     T=p/(rho*R);
     M=math.sqrt(u**2+v**2)/math.sqrt(abs(gamma*p/rho));
     Pt=p+0.5*rho*(u**2+v**2);
     Tt=T*(1+((gamma-1)/2))*M**2;
     Cp=(p-p_inf)/(0.5*rho_inf*u_inf**2);
     S=p/(abs(rho)**gamma);

     # copy now in post_proc field for this solution point
     post_proc[index][0] = u;
     post_proc[index][1] = v;
     post_proc[index][2] = p;
     post_proc[index][3] = T;
     post_proc[index][4] = M;
     post_proc[index][5] = Pt;
     post_proc[index][6] = Tt;
     post_proc[index][7] = Cp;
     post_proc[index][8] = S;


  mesh.write_mesh( file=cf.URI('file:euler-cylinder-2d-'+str(euler.time.current_time)+'.msh'),
                   fields=[solution.uri(),post_proc.uri(),euler.fields.dt.uri()])
  
#### Configure timestepping
#solver.Time.end_time = 3000000
#solver.Time.time_step = 1
#solver.TimeStepping.time_accurate = True          # time accurate for initial stability
##solver.iterative_solver = 'cf3.sdm.ExplicitRungeKuttaLowStorage2'
#solver.iterative_solver = 'cf3.sdm.lusgs.LUSGS'
#solver.TimeStepping.IterativeSolver.options.system = 'cf3.sdm.implicit.BDF2'
#solver.TimeStepping.cfl = 'min(4,0.05*(i+1.))'     # increasing cfl number
##solver.TimeStepping.cfl = '0.1'     # increasing cfl number
#solver.TimeStepping.max_iteration = 220           # limit the number of iterations (default = no limit)
#solver.TimeStepping.IterativeSolver.print_iteration_summary = True
#solver.TimeStepping.IterativeSolver.max_sweeps = 100
#solver.TimeStepping.IterativeSolver.convergence_level = 1e-7
#solver.TimeStepping.IterativeSolver.children.System.children.ComputeCellJacobian.options.reference_solution = [rho_inf,rho_inf*u_inf,rho_inf*u_inf,rhoE_inf]

##solver.TimeStepping.IterativeSolver.nb_stages = 3 # Runge Kutta number of stages


## Tecplot
##########
## Tecplot cannot write high-order meshes. A finer P1 mesh is loaded,
## and fields are interpolated to the P1-mesh. The mesh is finer to visualize
## the high-order solution better.

## Generate visualization mesh
#visualization_mesh = domain.load_mesh(file = cf.URI('../../../resources/cylinder-quad-p1-128x32.msh'), name = 'visualization_mesh');

## Interpolate fields using solution polynomial
#visualization_mesh.access_component('geometry').create_field(name='solution',  variables='rho[1],rhoU[2],rhoE[1]')
#visualization_mesh.access_component('geometry').create_field(name='wave_speed',variables='ws[1]')
#visualization_mesh.access_component('geometry').create_field(name='post_proc', variables='U[vec],p[1],T[1],M[1],Pt[1],Tt[1],Cp[1],S[1]')

#interpolator = model.get_child('tools').create_component('interpolator','cf3.mesh.actions.Interpolate')
#interpolator.interpolate(source=mesh.access_component("solution_space/solution").uri(),
#												 target=visualization_mesh.access_component("geometry/solution").uri())
#interpolator.interpolate(source=mesh.access_component("solution_space/wave_speed").uri(),
#													target=visualization_mesh.access_component("geometry/wave_speed").uri())
#interpolator.interpolate(source=mesh.access_component("solution_space/post_proc").uri(),
#												 target=visualization_mesh.access_component("geometry/post_proc").uri())

#fields = [
#visualization_mesh.access_component('geometry/solution').uri(),
#visualization_mesh.access_component('geometry/wave_speed').uri(),
#visualization_mesh.access_component('geometry/post_proc').uri()
#]

## Write visualization mesh
#tec_writer = model.get_child('tools').create_component('writer','cf3.mesh.tecplot.Writer')
#tec_writer.options().set('mesh',visualization_mesh)
#tec_writer.options().set('fields',fields)
#tec_writer.options().set('cell_centred',True)
#tec_writer.options().set('file',cf.URI('file:sdm_output.plt'))
#tec_writer.execute()
