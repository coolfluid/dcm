import coolfluid as cf

### Create new model specialized for SD
model   = cf.root.create_component('accousticpulse_2d','cf3.sdm.Model');

### Load the mesh
mesh = model.domain.load_mesh(file = cf.URI('../../../resources/circle-quad-p1-32.msh'), name = 'circle');

### Add the Partial Differential Equations to solve
lineuler = model.add_pde(name='lineuler',type='cf3.sdm.lineuler.LinEulerUniform2D',order=4)
lineuler.gamma = 1.4
lineuler.U0 = [0,0]
lineuler.rho0 = 1
lineuler.p0 = 1

### Add BC
lineuler.add_bc( name='farfield',
                 type='cf3.sdm.lineuler.BCFarfield2D',
                 regions=[ mesh.topology.boundary ] )

### Initialize the solution
model.tools.init_field.init_field(
  field=lineuler.solution,
  functions=[
    '0.001*exp( -( (x)^2 + (y)^2 )/(0.05)^2 )',
    '0',
    '0',
    '1.4 * 0.001*exp( -( (x)^2 + (y)^2 )/(0.05)^2 )' ] )

### Create the Solver for the Partial Differential Equations
solver = model.add_solver(pde=lineuler)
solver.children.time_step.cfl = 0.2
solver.children.scheme.nb_stages = 3

### Time Stepping
model.time_stepping.end_time = 0.3
model.time_stepping.execute()


#######################################
# POSTPROC to check accuracy
#######################################
exact_solution = lineuler.fields.create_field(name='exact_solution',variables='rho_ex[s],U_ex[v],p_ex[s]')
init_acousticpulse = model.tools.create_component('init_acousticpulse','cf3.sdm.lineuler.InitAcousticPulse')
init_acousticpulse.field = exact_solution
init_acousticpulse.time = 0.3
init_acousticpulse.execute()

solution = lineuler.fields.solution

difference = lineuler.fields.create_field(name='difference',variables='drho[s],dU[v],dp[s]')
for i in range(len(difference)) :
    difference[i][0] = exact_solution[i][0] - solution[i][0]
    difference[i][1] = exact_solution[i][1] - solution[i][1]/lineuler.rho0
    difference[i][2] = exact_solution[i][2] - solution[i][2]/lineuler.rho0
    difference[i][3] = exact_solution[i][3] - solution[i][3]

compute_norm = model.tools.create_component('compute_norm','cf3.sdm.ComputeLNorm')
compute_norm.field = difference
compute_norm.order = 2
compute_norm.execute()
print "norm = ",compute_norm.properties()['norm']





########################
# OUTPUT
########################

fields = [
  lineuler.fields.solution.uri(),
  lineuler.fields.ws.uri(),
  lineuler.fields.rhs.uri(),
  lineuler.fields.exact_solution.uri(),
  lineuler.fields.difference.uri()
]

mesh.write_mesh(file=cf.URI('file:lineuler-acousticpulse-2d.plt'),fields=fields)
mesh.write_mesh(file=cf.URI('file:lineuler-acousticpulse-2d.msh'),fields=fields)

