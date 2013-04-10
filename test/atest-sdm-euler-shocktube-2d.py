import coolfluid as cf

cf.env.log_level = 3

### Create new model
model = cf.root.create_component('shocktube_2d','cf3.sdm.Model');

### Generate mesh
mesh = model.domain.create_component('mesh','cf3.mesh.Mesh')
mesh_generator = model.tools.create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")
mesh_generator.options().set("mesh",mesh.uri())
mesh_generator.options().set("nb_cells",[20,20])
mesh_generator.options().set("lengths",[10,10])
mesh_generator.options().set("offsets",[-5,-5])
mesh_generator.execute()
load_balance = mesh_generator.create_component("load_balancer","cf3.mesh.actions.LoadBalance")
load_balance.options().set("mesh",mesh)
load_balance.execute()

### Add PDE
euler = model.add_pde(name='euler',
                      type='cf3.sdm.equations.navierstokes.NavierStokes2D',
                      shape_function='cf3.sdm.core.LegendreGaussEndP2')
euler.gamma = 1.4
euler.R = 287.05

euler.add_bc(name='mirror',type='cf3.sdm.equations.euler.BCMirror2D',regions=
[ mesh.topology.left, mesh.topology.right, mesh.topology.top, mesh.topology.bottom ])

### Initial solution
model.tools.init_field.init_field( field=euler.solution, functions = 
[ 'r_L:=4.696; r_R:=1.408; if(x<=0 & y<=0,r_L,r_R)',
  '0.','0.',
  'p_L:=404400; p_R:=101100; g:=1.4; if(x<=0 & y<=0,p_L/(g-1),p_R/(g-1))'
] )

### Solve
solver = model.add_solver(pde=euler,solver='cf3.sdm.solver.erk.MidPoint')
solver.children.time_step_computer.cfl=0.01
solver.solve_time_step(0.008)

mesh.write_mesh(file=cf.URI('file:euler-shocktube-2d.msh'), fields=[euler.solution.uri()])


