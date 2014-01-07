from util import log, read, write, symlink, mkdir
from coolfluid import URI

def save(simulation):
    log( 'Saving ',simulation.name+'_iter'+str(simulation.pde.time.iteration) )

    case_dir = simulation.name+'.case'
    case_iter_dir = case_dir+'/iter'+str(simulation.pde.time.iteration).zfill(8)+'.case'
    mkdir(case_dir)
    mkdir(case_iter_dir)
    simulation.write_mesh( case_iter_dir+'/mesh.cf3mesh' )
    write( simulation, case_iter_dir+'/case.p' )
    simulation.solver.children.history.write(URI(case_iter_dir+'/history.tsv'))
    symlink( 'iter'+str(simulation.pde.time.iteration).zfill(8)+'.case', case_dir+'/current.case' )
    
def restart(case_dir,iteration='current'):
    
    if isinstance(iteration,str):
        case_iter_dir=case_dir+'/'+iteration+'.case'
    else:
        case_iter_dir=case_dir+'/iter'+str(iteration).zfill(8)+'.case'

    simulation = read(case_iter_dir+'/case.p')
    log('Restarting ',case_iter_dir)
    simulation.read_mesh(case_iter_dir+'/mesh.cf3mesh')
    simulation.create_space_discretization( order = simulation.order, riemann_solver = simulation.riemann_solver )
    simulation.pde.time.current_time = simulation.time
    simulation.pde.time.iteration = simulation.iteration
    bcs = simulation.bcs
    simulation.bcs = []
    for bc in bcs:
        simulation.add_bc(str(bc[0]),str(bc[1]),bc[2],**bc[3])
    simulation.set_time_discretization( type = simulation.solver_type )
    simulation.set_cfl(str(simulation.max_cfl)) # set the current cfl number
    simulation.set_cfl(str(simulation.cfl)) # set the cfl function
    simulation.set_time_accurate(simulation.time_accurate)
    simulation.solver.children.history.read(URI(case_iter_dir+'/history.tsv'))
    return simulation
