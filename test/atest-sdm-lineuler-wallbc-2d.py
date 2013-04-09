import coolfluid as cf

### Create new model specialized for SD
model   = cf.root.create_component('accousticpulse_2d','cf3.sdm.Model');

### Load the mesh
mesh = model.domain.create_component('mesh','cf3.mesh.Mesh')
mesh_gen = model.tools.create_component('mesh_gen','cf3.mesh.SimpleMeshGenerator')
mesh_gen.mesh = mesh.uri()
mesh_gen.lengths = [200,200]
mesh_gen.offsets = [-100,0]
mesh_gen.nb_cells = [30,30]
mesh_gen.execute()

repartitioner=mesh.create_component('repartitioner','cf3.mesh.actions.LoadBalance')
repartitioner.mesh = mesh
repartitioner.execute()

### Add the Partial Differential Equations to solve
lineuler = model.add_pde(name='lineuler',type='cf3.sdm.equations.lineuler.LinEulerUniform2D',
    shape_function='cf3.sdm.core.LegendreGaussLobattoP2')
lineuler.gamma = 1.
lineuler.U0 = [0.5,0]
lineuler.rho0 = 1
lineuler.p0 = 1

### Add BC
lineuler.add_bc( name='farfield', type='cf3.sdm.equations.lineuler.BCFarfield2D',
                 regions=[ mesh.topology.left, mesh.topology.top ] )
lineuler.add_bc( name='mirror', type='cf3.sdm.equations.lineuler.BCMirror2D',
                 regions=[ mesh.topology.bottom ] )
lineuler.add_bc( name='outlet', type='cf3.sdm.equations.lineuler.BCExtrapolation2D',
                 regions=[ mesh.topology.right ] )

### Initialize the solution
model.tools.init_field.init_field(
  field=lineuler.solution,
  functions=[
    'exp( -log(2)*( (x)^2 + (y-25)^2 )/(25) )',
    '0',
    '0',
    '1.^2 * exp( -log(2)*( (x)^2 + (y-25)^2 )/(25) )' ] )

### Create the Solver for the Partial Differential Equations
solver = model.add_solver(pde=lineuler)
solver.children.time_step_computer.cfl = 0.3
solver.options.order = 3

### Time Stepping
model.time_stepping.end_time = 100
model.time_stepping.execute()

#######################################
# POST-PROCESSING
#######################################

compute_char = model.create_component('compute_characteristics','cf3.sdm.equations.lineuler.ComputeCharacteristicVariablesUniform2D')
compute_char.normal = [0.,-1.]
compute_char.field = lineuler.solution
compute_char.c0 = 1.
compute_char.execute()

########################
# OUTPUT
########################

fields = [
lineuler.fields.solution.uri(),
lineuler.fields.char.uri()
]

vis_mesh = model.domain.create_component('vis_mesh','cf3.mesh.Mesh')
mesh_gen.mesh = vis_mesh.uri()
mesh_gen.nb_cells = [100,100]
mesh_gen.execute()
vis_solution = vis_mesh.geometry.create_field(name='solution',variables='rho,rho0u[vector],p')

interpolator = model.tools.create_component('interpolator','cf3.mesh.ShapeFunctionInterpolator')
interpolator.interpolate(source=lineuler.solution.uri(),target=vis_solution.uri())

vis_mesh.write_mesh(file=cf.URI('file:wallbc.plt'), fields=[vis_solution.uri()])
vis_mesh.write_mesh(file=cf.URI('file:wallbc.msh'), fields=[vis_solution.uri()])
