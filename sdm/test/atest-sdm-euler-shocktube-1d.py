import coolfluid as cf

cf.env.log_level = 3

### Create new model
model = cf.root.create_component('shocktube_1d','cf3.dcm.Model');

### Generate mesh
mesh = model.domain.create_component('mesh','cf3.mesh.Mesh')
mesh_generator = model.tools.create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")
mesh_generator.mesh = mesh.uri()
mesh_generator.nb_cells = [70]
mesh_generator.lengths = [10]
mesh_generator.offsets = [-5]
mesh_generator.execute()
load_balance = mesh_generator.create_component("load_balancer","cf3.mesh.actions.LoadBalance")
load_balance.options().set("mesh",mesh)
load_balance.execute()
model.build_faces();

### Add PDE
euler = model.add_pde(name='euler',type='cf3.dcm.equations.euler.Euler1D',shape_function='cf3.dcm.core.LegendreGaussEndP2')
euler.add_term(name='rhs',type='cf3.sdm.br2.euler_RightHandSide1D')
euler.gamma = 1.4
euler.R = 287.05

euler.add_bc(name='reflect',type='cf3.dcm.equations.euler.BCMirror1D',regions=[mesh.topology.xneg,mesh.topology.xpos])

### Initial solution
model.tools.init_field.init_field( field=euler.solution, functions = 
[ 'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,r_L,r_R)',
  'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,r_L*u_L,r_R*u_R)',
  'r_L:=4.696; r_R:=1.408; u_L:=0; u_R:=0; p_L:=404400; p_R:=101100; g:=1.4; if(x<=0,p_L/(g-1)+0.5*r_L*u_L*u_L,p_R/(g-1)+0.5*r_R*u_R*u_R)'
] )

### Solve
solver = model.add_solver(name='solver',pde=euler,solver='cf3.dcm.solver.erk.MidPoint')
solver.children.time_step_computer.cfl=0.3
solver.solve_time_step(0.008)

mesh.write_mesh(file=cf.URI('file:euler-shocktube-1d.plt'), fields=[euler.solution.uri()])

