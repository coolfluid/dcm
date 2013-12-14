import coolfluid as cf
import math
### Create new model specialized for SD
model   = cf.root.create_component('accousticpulse_2d','cf3.dcm.Model');

### Load the mesh
mesh = model.domain.load_mesh(file = cf.URI('../../../resources/square100-quad-p2-50x50.msh'), name = 'square')
model.build_faces();

### Add the Partial Differential Equations to solve
lineuler = model.add_pde(name='lineuler',type='cf3.dcm.equations.lineuler.LinEuler2D',
    shape_function='cf3.dcm.core.LegendreGaussEndP3')
lineuler.gamma = 1.
U0 = [0.5,0]
rho0 = 1
p0 = 1

lineuler.add_term( name='rhs', type='cf3.sdm.br2.lineuler_RightHandSide2D' )

### Add BC
lineuler.add_bc( name='farfield',
                 type='cf3.dcm.equations.lineuler.BCFarfield2D',
                 regions=[ mesh.topology.left, mesh.topology.bottom, mesh.topology.top ] )
lineuler.add_bc( name='outlet',
                 type='cf3.dcm.equations.lineuler.BCExtrapolation2D',
                 regions=[ mesh.topology.right ] )

### Initialize the solution
model.tools.init_field.init_field(
  field=lineuler.solution,
  functions=[
    ' exp( -log(2.)*((x)^2+y^2)/9. ) + 0.1*exp( -log(2.)*((x-67.)^2 + y^2)/25. )',
    ' 1.*0.04*y      *exp( -log(2.)*((x-67.)^2+y^2)/25. )',
    '-1.*0.04*(x-67.)*exp( -log(2.)*((x-67.)^2+y^2)/25. )',
    '1.* exp( -log(2.)*((x)^2+y^2)/9. )' ] )
model.tools.init_field.init_field(
  field=lineuler.background,
  functions=[ str(rho0), str(U0[0]), str(U0[1]), str(p0) ] )

model.tools.init_field.init_field(
  field=lineuler.bdry_background,
  functions=[ str(rho0), str(U0[0]), str(U0[1]), str(p0) ] )

### Create the Solver for the Partial Differential Equations
solver = model.add_solver(pde=lineuler,name='optim_erk',solver='cf3.sdm.solver.optim_erkls.ERK_18_4')
solver.children.time_step_computer.cfl = 1.5*1.17418695241

### Time Stepping
model.time_stepping.end_time = 1 #90
model.time_stepping.time_step = 1 #10

while not model.time_stepping.properties.finished :
    model.time_stepping.do_step()
    mesh.write_mesh(file=cf.URI('solution'+str(model.time_stepping.step)+'.msh'),fields=[lineuler.solution.uri()])

## function describing entropy and vortex without acoustic pulse
#entropy_vortex = [
# '0.1*exp( -log(2.)*((x-67)^2 + y^2)/25. )',
# ' 0.04*y      *exp( -log(2.)*((x-67)^2+y^2)/25. )',
# '-0.04*(x-67)*exp( -log(2.)*((x-67)^2+y^2)/25. )',
# '0'
#]

## function describing acoustic pulse only
#acoustic = [
# 'exp( -log(2.)*((x)^2+y^2)/9. )',
# ' 0',
# '-0',
# 'exp( -log(2.)*((x)^2+y^2)/9. )'
#]


#######################################
# POST-PROCESSING
#######################################

compute_char = model.tools.create_component('compute_characteristics','cf3.dcm.equations.lineuler.ComputeCharacteristicVariablesUniform2D')
compute_char.options().set('normal',[1.,0.])
compute_char.options().set('field',lineuler.solution)
compute_char.options().set('c0',math.sqrt(lineuler.gamma*p0/rho0))
compute_char.execute()

########################
# OUTPUT
########################

fields = [
lineuler.fields.solution.uri(),
lineuler.fields.char.uri(),
lineuler.fields.gradn_char.uri(),
]

mesh.write_mesh(file=cf.URI('file:lineuler-acousticvorticity-2d.msh'),fields=fields)


# Tecplot
#########
# Tecplot cannot write high-order meshes. A finer P1 mesh is generated,
# and fields are interpolated to the P1-mesh. The mesh is finer to visualize
# the high-order solution better.

mesh_generator = model.tools.create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")

# Generate visualization mesh
visualization_mesh = model.domain.create_component('visualization_mesh','cf3.mesh.Mesh')
mesh_generator.options().set("mesh",visualization_mesh.uri())
mesh_generator.options().set("nb_cells",[400,400])
mesh_generator.options().set("lengths",[200,200])
mesh_generator.options().set("offsets",[-100,-100])
mesh_generator.execute()

# Interpolate fields using solution polynomial
visualization_mesh.geometry.create_field(name='solution',       variables='rho[1],rho0U[2],p[1]')
#visualization_mesh.get_child('geometry').create_field(name='char', variables='S[1],Shear[1],Aplus[1],Amin[1],A[1],omega[1]')

interpolator = model.tools.create_component('interpolator','cf3.mesh.ShapeFunctionInterpolator')
interpolator.interpolate(source=lineuler.fields.solution.uri(),
                         target=visualization_mesh.geometry.solution.uri())
#interpolator.interpolate(source=mesh.access_component("solution_space/char").uri(),
#												 target=visualization_mesh.access_component("geometry/char").uri())

fields = [
visualization_mesh.geometry.solution.uri(),
#visualization_mesh.access_component('geometry/char').uri()
]

# Write visualization mesh
visualization_mesh.write_mesh(file=cf.URI('file:lineuler-acousticvorticity-2d.plt'),fields=fields)

#####################
# Probe line y=0
#####################

# Generate 1D line mesh, for now only y=0 can be probed as the line has 1D coordinates only
probe_mesh = model.domain.create_component('probe_mesh','cf3.mesh.Mesh')
mesh_generator.options().set("mesh",probe_mesh.uri())
mesh_generator.options().set("nb_cells",[1000])
mesh_generator.options().set("lengths",[200])
mesh_generator.options().set("offsets",[-100])
mesh_generator.execute()

# Interpolate fields
probe_mesh.get_child('geometry').create_field(name='solution', variables='rho[1],rho0U[2],p[1]')
#probe_mesh.get_child('geometry').create_field(name='char',     variables='S[1],Shear[1],Aplus[1],Amin[1],A[1],omega[1]')

interpolator.interpolate(source=lineuler.fields.solution.uri(),
                         target=probe_mesh.geometry.solution.uri())
#interpolator.interpolate(source=mesh.access_component("solution_space/char").uri(),
#												 target=probe_mesh.access_component("geometry/char").uri())

fields = [
probe_mesh.geometry.solution.uri(),
#probe_mesh.access_component('geometry/char').uri()
]

# Write probe mesh
probe_mesh.write_mesh(file=cf.URI('file:probe_liney0.plt'),fields=fields)

