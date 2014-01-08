#! /usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser( prog="mpirun -np TASKS output.py", description=
"""
Read DCM case, post-process, and output mesh.
Must be invoked through mpirun, and restarted with same number of tasks as before
""")
parser.add_argument('case',  type=str, help='DCM case')
parser.add_argument('-f', '--fields', type=str, nargs='+', help='space-separated list of field names (default=solution)')
parser.add_argument('-o', '--output', type=str, help='Output mesh (default=out.msh)', default='out.msh')
parser.add_argument('--iteration',    type=int, help='iteration to process (default=latest)')
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

if args.fields:
    log("Computing fields")
    for field in args.fields:
        log("  "+field)
    sim.postprocessing()

sim.write_mesh(args.output,args.fields)
