#! /usr/bin/env python

"""
Case:
-----
    2D Navier Stokes simulation around a square cylinder
    Reynolds number: 10000
    Mach: 0.2
    
Mesh:
-----

Create following mesh using gmsh
>>> gmsh square-cylinder.geo

          farfield
 __________________________________
|                                  |
| farfield                         |
|      ___    o     o              |
|     |   |  o  o  o  o            |
|     |___| o     o    o           |
|      cylinder                    | outlet
|                                  |
|__________________________________|
            farfield

"""


from pydcm import NavierStokes, log
from pydcm import save, restart

sim = NavierStokes( 'square_cylinder', dimension = 2 )

# Read mesh
sim.read_mesh('square-cylinder.msh')

sim.define( gamma = 1.4 )
sim.define( h = 1., Re = 1e4 ,  Mref = 0.2 , R = 287.05, Pr = 0.72 )
sim.define( pref = 1.e5, rhoref = 1.)
sim.define( cref = 'sqrt(gamma*pref/rhoref)' )
sim.define( Uref = 'Mref*cref' )
sim.define( nu = 'Uref*h/Re' )
sim.define( mu = 'rhoref*nu' )
sim.define( Cp = 'gamma*R/(gamma-1)' )
sim.define( kappa = 'Cp*mu/Pr' )
sim.define( Eref = 'pref/(gamma-1)/rhoref + 0.5*(Uref**2)' )
sim.define( Tref = 'pref / (rhoref * R)' )

sim.print_constants()

# Choose the space discretization order
sim.create_space_discretization( order = 1 )

# Add boundary conditions
sim.add_bc( 'farfield',  'Farfield', ['farfield'], rho='rhoref', U=['Uref','0'], p='pref' )
sim.add_bc( 'outflow', 'PressureOutlet', ['outlet'], p='pref' )
sim.add_bc( 'wall',   'AdiabaticWall', ['cylinder'] )

sim.set_time_discretization( 'LUSGS_BDF2' )
sim.set_time_accurate( True )
sim.set_cfl('1.')
sim.set_cfl('min(10, cfl*1.01)')

sim.init_solution(
    'rhoref',
    'rhoref*Uref',
    '0',
    'rhoref*Eref' )

while True:
    sim.propagate(iterations=1000)
    save(sim)
