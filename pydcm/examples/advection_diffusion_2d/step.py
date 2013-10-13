from pydcm import AdvectionDiffusion, log

sim = AdvectionDiffusion( 'step', dimension = 2 )

# Choose either to generate a mesh, or another mesh
# sim.create_simple_mesh( lengths=[110.,20.], nb_cells=[20,4], offsets=[-10.,0.] )
sim.read_mesh('meshes/rectangle.msh')

# Set advection velocity and diffusion coefficient
sim.define( a = [1.,0.] ,  mu = .25 )

# Choose the space discretization order
sim.create_space_discretization( order = 5 )

# Add boundary conditions
sim.add_bc( 'inlet',  'Dirichlet',     ['left'],         Q = 5.    )
sim.add_bc( 'sides',  'Neumann',       ['bottom','top'], dQdn = 0. )
sim.add_bc( 'outlet', 'Extrapolation', ['right']                   )

# Choose the time discretization method
sim.set_time_discretization( 'ERK_18_4' )
sim.set_cfl(1.2)

sim.init_solution( 'if(x<=0., 5., -5.)' )

sim.propagate( time_step=50. )

sim.write_mesh('output/step_50.msh')
