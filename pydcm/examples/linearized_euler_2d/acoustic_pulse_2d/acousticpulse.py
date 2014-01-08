#! /usr/bin/env python

"""
Case:
-----
    2D Linearized Euler simulation of Acoustic pulse
    
Create following mesh using gmsh
>>> gmsh circle-quad.geo
"""
from pydcm import LinearizedEuler, log

sim = LinearizedEuler( 'acousticpulse', dimension = 2 )

# Choose either to generate a mesh, or read another mesh
sim.create_simple_mesh( lengths=[1.,1.], nb_cells=[31,31], offsets=[-0.5,-0.5] )
bdry = ['left','bottom','right','top']

# sim.read_mesh('meshes/circle-quad.msh')
# bdry = ['boundary']

sim.define( gamma = 1.4 )
sim.define( rho0 = 1., U0 = [0.,0.] ,  p0 = 1. )
sim.define( c0 = 'sqrt(p0*gamma/rho0)' )

# Choose the space discretization order
sim.create_space_discretization( order = 3 )

# Add boundary conditions

sim.add_bc( 'farfield',  'Farfield', bdry )

# Choose the time discretization method
sim.set_time_discretization( 'MidPoint' )
sim.set_cfl(0.5)

sim.init_solution( 
    '0.001*exp( -( (x)^2 + (y)^2 )/(0.05)^2 )',
    '0',
    '0',
    'c0^2 * 0.001*exp( -( (x)^2 + (y)^2 )/(0.05)^2 )' )
    
sim.init_background( 'rho0', 'U0[x]', 'U0[y]', 'p0' )

sim.fields.append( sim.pde.fields.background.uri() )
sim.fields.append( sim.pde.fields.background_gradient.uri() )

sim.write_mesh( 'output/acousticpulse_00.msh' )

sim.propagate( time_step = 0.3 )

sim.write_mesh( 'output/acousticpulse_03.msh' )
