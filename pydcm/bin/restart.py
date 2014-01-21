#! /usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser( prog="mpirun -np TASKS restart.py", description=
"""
Restart DCM case, with several changes in configuration.
Must be invoked through mpirun, and restarted with same number of tasks as before
""")
parser.add_argument('case',  type=str, help='DCM case')
parser.add_argument('--iteration', type=int, help='Iteration to restart from (default=latest)')
parser.add_argument('--save', type=int, help='Save solution every SAVE iterations (default=1000)', default=1000)
parser.add_argument('-o', '--order', type=int, help='Change the order of simulation')
parser.add_argument('--cfl', type=int, help='Change the CFL number gradually to this number')
args = parser.parse_args()


from pydcm import restart, save, log

if args.iteration:
    iteration = args.iteration
else:
    iteration = 'current'

sim = restart(args.case,iteration=iteration)
log("Loaded simulation "+sim.name)
log("  iteration: "+str(sim.iteration))
log("  time:      "+str(sim.time))

sim.print_constants()

if args.order:
    sim.create_space_discretization( order=args.order )

if args.cfl:
    sim.set_cfl('if(cfl<=0,1,min('+str(args.cfl)+',cfl*1.01))')

sim.solver.convergence_level = 1e-6

from time import gmtime, strftime
while True:
    log("Started Computations at "+strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    sim.propagate( iterations=args.save )
    save(sim)