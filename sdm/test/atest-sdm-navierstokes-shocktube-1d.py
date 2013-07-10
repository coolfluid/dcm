import coolfluid as cf

cf.env.log_level = 3

### Create new model
model = cf.root.create_component('shocktube_1d','cf3.sdm.Model');

### Generate mesh
mesh = model.domain.create_component('mesh','cf3.mesh.Mesh')
mesh_generator = model.tools.create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")
mesh_generator.options().set("mesh",mesh.uri())
mesh_generator.options().set("nb_cells",[70])
mesh_generator.options().set("lengths",[10])
mesh_generator.options().set("offsets",[-5])
mesh_generator.execute()
load_balance = mesh_generator.create_component("load_balancer","cf3.mesh.actions.LoadBalance")
load_balance.options().set("mesh",mesh)
load_balance.execute()

### Add PDE
# navierstokes = model.add_pde(name='navierstokes',type='cf3.sdm.equations.euler.Euler1D',order=1)
navierstokes = model.add_pde(name='navierstokes',type='cf3.sdm.equations.navierstokes.NavierStokes1D',order=1)
navierstokes.gamma = 1.4
navierstokes.R = 287.05
navierstokes.mu = 25

### Initial solution
model.tools.init_field.init_field( field=navierstokes.solution, functions = 
[ 'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,r_L,r_R)',
  'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,r_L*u_L,r_R*u_R)',
  'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,p_L/(g-1)+0.5*r_L*u_L*u_L,p_R/(g-1)+0.5*r_R*u_R*u_R)'
] )

### Solve
solver = model.add_solver(pde=navierstokes,solver='cf3.sdm.solver.erk.ForwardEuler')
solver.children.time_step_computer.cfl=0.3
solver.solve_time_step(0.008)

mesh.write_mesh(file=cf.URI('file:navierstokes-shocktube-1d.plt'), fields=[navierstokes.solution.uri()])

