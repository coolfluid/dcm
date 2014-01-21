from pydcm.util import Definitions, log
import coolfluid as cf
import os


class DCM(object):
    """
    Spectral Difference Method simulation

    Implements basic model structure.
    Model subroutines should be called in a certain order:
    - define( var = expression )  For as many variables that are used
    - read_mesh(filename)
    - create_space_discretization(order,riemann_solver)
    - add_bc(name, type, regions, extra_bc_config )
    - set_time_discretization(type)
    - set_cfl( expression(i,t,cfl) )
    - set_time_accurate( true or false )
    - propagate( iterations=nb_iterations )
    """
    

    def __init__(self,name='sdm',dimension=0):
        self.name = name
        self.dimension = dimension
        self.model = cf.root.create_component(name,'cf3.dcm.Model');
        self.shape_function = 'cf3.dcm.core.LegendreGaussEnd'
        self.fields = []
        self.__time_accurate = True
        self.__defs = Definitions()
        self.order = 0
        self.solver_type = 'LUSGS_BDF2'
        self.riemann_solver = 'Roe'
        self.cfl = 1.
        self.bcs = []
        self.pde = None
       
    @property 
    def time_accurate(self):
        return self.__time_accurate

    @property
    def time(self):
        return self.pde.time.current_time

    def read_mesh(self,filename):
        if self.dimension != 0:
            self.model.domain.options.dimension = self.dimension

        if ( not os.path.exists(filename) ):
            raise RuntimeError, "Could not load mesh because "+filename+" doesn't exist"

        if (os.path.splitext(filename)[1] == '.cf3mesh' ):
            reader = self.model.tools.create_component('cf3mesh_reader', 'cf3.mesh.cf3mesh.Reader')
            self.mesh = reader.read( domain=self.model.domain, file=cf.URI(filename) )            
        else:
            self.mesh = self.model.domain.load_mesh( file=cf.URI(filename), name='Mesh' );
            self.model.build_faces()
            
        self.dimension = self.mesh.properties.dimension
        
    def create_simple_mesh(self,lengths,nb_cells,offsets=[]):
        if self.dimension != 0:
            self.model.domain.options.dimension = self.dimension
        self.mesh = self.model.domain.create_component('Mesh','cf3.mesh.Mesh')
        mesh_generator = self.model.tools.create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")
        mesh_generator.mesh = self.mesh.uri()
        mesh_generator.nb_cells = nb_cells
        mesh_generator.lengths = lengths
        if offsets: mesh_generator.offsets = offsets
        mesh_generator.execute()
        load_balance = mesh_generator.create_component("load_balancer","cf3.mesh.actions.LoadBalance")
        load_balance.options().set("mesh",self.mesh)
        load_balance.execute()
        self.model.build_faces()
    

    def write_mesh(self,filename,fields=[]):
        if not fields:
            if self.fields:
                self.mesh.write_mesh(file=cf.URI(filename), fields=self.fields)
            else:
                self.mesh.write_mesh(file=cf.URI(filename))
        else:
            fields_uri = []
            for field in fields:
                if ( self.pde.fields.get_child(field) ):
                    fields_uri.append( self.pde.fields.get_child(field).uri() )
                else:
                    raise RuntimeError, "Could not find field "+field
            self.mesh.write_mesh(file=cf.URI(filename), fields=fields_uri)
            
        
    def create_space_discretization(self,order,riemann_solver=''):
        if (self.pde is not None):
            if (self.order == order):
                return
            prev_pde = self.pde
            prev_order = self.order
            prev_fields = self.fields
            prev_bcs = self.bcs
            self.fields = []
            self.bcs = []
            
            self.order=order
            self.riemann_solver=riemann_solver
            self.pde = self.model.add_pde( 
                name='pde_P'+str(self.order-1),
                type=self.pde_type+str(self.dimension)+'D',
                shape_function=self.shape_function+'P'+str(self.order-1) )
            if riemann_solver: self.pde.riemann_solver = riemann_solver
            if not (str(self.pde.solution) in [ str(f) for f in self.fields] ):
                self.fields.append(self.pde.solution.uri())
            self.add_default_terms()

            
            self.pde.time.current_time = prev_pde.time.current_time
            self.pde.time.iteration = prev_pde.time.iteration

            for bc in prev_bcs:
                self.add_bc(str(bc[0]),str(bc[1]),bc[2],**bc[3])
            self.set_time_discretization( type = self.solver_type )
            self.set_cfl(str(self.cfl))
            self.set_time_accurate(self.time_accurate)
            
            self.interpolate_solution(prev_pde.solution,self.pde.solution)
            
            prev_pde.delete_component()
            
        else:            
            self.order=order
            self.riemann_solver=riemann_solver
            self.pde = self.model.add_pde( 
                name='pde_P'+str(self.order-1),
                type=self.pde_type+str(self.dimension)+'D',
                shape_function=self.shape_function+'P'+str(self.order-1) )
            if riemann_solver: self.pde.riemann_solver = riemann_solver
            if not (str(self.pde.solution) in [ str(f) for f in self.fields] ):
                self.fields.append(self.pde.solution.uri())
            self.add_default_terms()
            
    def set_time_discretization(self,type='LUSGS_BDF2'):
        self.solver_type = type
        if (type == 'LUSGS_BDF1'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.sdmx.solver.lusgs.BDF1')
            self.solver.jacobian.reference_solution = self.ref_solution()
            self.solver.convergence_level = 1e-8
            self.solver.recompute_lhs_frequency = 5
            self.solver.min_sweeps=1
            self.solver.max_sweeps=20
            self.solver.print_sweep_summary=False
            self.solver.children.time_step_computer.cfl = 'if(i<100,max( 0.1 , min(cfl*1.1, 30) ), max(0.1,min(cfl*1.01, 500)))'
            self.solver.adaptive_time_stepping = False
        elif (type == 'LUSGS_BDF2'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.sdmx.solver.lusgs.BDF2')
            self.solver.jacobian.reference_solution = self.ref_solution()
            self.solver.convergence_level = 1e-8
            self.solver.recompute_lhs_frequency = 5
            self.solver.min_sweeps=1
            self.solver.max_sweeps=20
            self.solver.print_sweep_summary=False
            self.solver.children.time_step_computer.cfl = 'if(i<100,max( 0.1 , min(cfl*1.1, 30) ), max(0.1,min(cfl*1.01, 500)))'
            self.solver.adaptive_time_stepping = False
        elif (type == 'ERK_17_3'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.sdm.solver.optim_erkls.ERK_17_3')
            self.set_cfl(1.82207767386)
        elif (type == 'ERK_18_4'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.sdm.solver.optim_erkls.ERK_18_4')
            self.set_cfl(1.17418695241)
        elif (type == 'ERK_20_5'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.sdm.solver.optim_erkls.ERK_20_5')
            self.set_cfl(0.843896716833)
        elif (type == 'ForwardEuler'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.ForwardEuler')
        elif (type == 'MidPoint'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.MidPoint')
        elif (type == 'Heun2'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.Heun2')
        elif (type == 'Heun3'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.Heun3')
        elif (type == 'ClassicRK33'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.ClassicRK33')
        elif (type == 'ClassicRK44'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.ClassicRK44')
        elif (type == 'LowStorageRK3'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erkls.TwoSstar')
            self.solver.order = 3
        elif (type == 'LowStorageRK4'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erkls.TwoSstar')
            self.solver.order = 4
        elif (type == 'SSPRK54'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.SSPRK54')
        elif (type == 'RK65'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.RK65')
        elif (type == 'RKF65'):
            self.solver = self.model.add_solver(name=type,pde=self.pde,solver='cf3.dcm.solver.erk.RKF65')
        else:
            raise RuntimeError, 'Time discretization type '+type+' not recognised'
        self.solver.children.time_step_computer.time_accurate = self.__time_accurate
        self.solver.children.history.file = cf.URI(self.name+'_history.tsv')
        
        self.set_post_iteration()

    def set_post_iteration(self):
        pass

    def save_bc(self,name,type,regions,**keyword_args):
        self.bcs.append( [name,type,regions,keyword_args] )
    
    def add_bc(self, name, type, regions, **keyword_args ):
        self.save_bc(name,type,regions,**keyword_args)        
        region_comps = [ self.mesh.topology.access_component(str(reg)) for reg in regions ]                
        bc = self.pde.add_bc( name=name, type=type, regions=region_comps )
        return bc

    def set_cfl(self, expression):
        self.solver.children.time_step_computer.cfl = expression
        self.cfl = expression

    def set_time_accurate(self, time_accurate=True):
        self.__time_accurate = time_accurate
        self.solver.children.time_step_computer.time_accurate = time_accurate
        if not time_accurate:
            self.pde.time.end_time = 1.
        else:
            self.solver.options.max_iteration = 0

    def define(self,**kwargs):
        self.__defs.set(**kwargs)
        for (key,value) in kwargs.iteritems():
            setattr(self,key,value)
            
    def defined(self,var):
        return self.__defs.has(var)

    def value(self,var):
        return self.__defs.eval(var)
        
    def expression(self,expr=0):
        return self.__defs.fparser(expr)

    def init_solution(self,*functions):
        self.model.tools.init_field.init_field(
          field=self.pde.solution,
          functions= [self.expression(func) for func in functions] )

    def interpolate_solution(self,source,target):
        interpolator = self.model.tools.get_child('interpolator')
        if ( not interpolator ) :
            interpolator = self.model.tools.create_component('interpolator','cf3.mesh.SpaceInterpolator')        
        interpolator.interpolate( source=source.uri(), 
                                  target=target.uri() )
        
    def propagate(self, **keyword_args):
        if 'end_time' in keyword_args:
            self.pde.time.end_time = keyword_args['end_time']

        if self.time_accurate:
            if 'iterations' in keyword_args:
                self.solver.options.max_iteration = self.pde.time.iteration + keyword_args['iterations']
                if not 'time_step' in keyword_args or not 'end_time' in keyword_args:
                    self.pde.time.end_time = 1e300
            else:
                self.solver.options.max_iteration = int(1e9)
            if 'time_step' in keyword_args:
                self.solver.solve_time_step( time_step = keyword_args['time_step'] )
            else:
                self.solver.execute()
        else:
            if self.pde.time.end_time == 0:
                self.pde.time.end_time = 1
            self.solver.solve_iterations(iterations=keyword_args['iterations'] )
            
    def postprocessing(self):
        pass
            
    def __getstate__(self):
        state = {}
        state['name'] = self.name
        state['pde_type']  = self.pde_type
        state['solver_type'] = self.solver_type
        state['order']       = self.order
        state['riemann_solver'] = self.riemann_solver
        state['dimension']   = self.dimension
        state['iteration']   = self.pde.time.iteration
        state['time']        = self.pde.time.current_time
        state['definitions']  = self.__defs
        state['fields'] = []
        for field in self.fields:
            state['fields'].append(str(field))
        state['time_accurate'] = self.time_accurate
        state['max_cfl'] = self.solver.children.time_step_computer.max_cfl()
        state['cfl'] = self.cfl
        state['bcs'] = self.bcs
        return state
    
    def __setstate__(self,state):
        self.name = str(state['name'])
        self.dimension = state['dimension']
        self.__defs = state['definitions']
        self.pde = None
        self.pde_type = str(state['pde_type'])
        self.order = state['order']
        self.riemann_solver = str(state['riemann_solver'])
        self.solver_type = str(state['solver_type'])
        self.model = cf.root.create_component(self.name,'cf3.dcm.Model');
        self.shape_function = str('cf3.dcm.core.LegendreGaussEnd')
        self.fields = []
        for field in state['fields']:
            self.fields.append( cf.URI(str(field) ) )
        self.__time_accurate = state['time_accurate']
        self.max_cfl = state.get('max_cfl',1.)
        self.cfl = state['cfl']
        self.begin_time = state['time']
        self.iteration = state['iteration']
        self.bcs = state['bcs']
        
        
