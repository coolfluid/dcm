from pydcm.dcm import DCM
from pydcm.util import log
import math

class NavierStokes(DCM):
    
    def __init__(self,name='navierstokes',dimension=0):
        super(NavierStokes, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.navierstokes.NavierStokes'
        self.turb_stats_count = -1 # neg value means OFF
        self.bdry_layer_regions = ['.']
                        
    def ref_solution(self):
        rho  = self.value('rhoref')
        rhoU = rho*self.value('Uref')
        rhoE = rho*self.value('Eref')
        return { 1: [ rho, rhoU, rhoE ],
                 2: [ rho, rhoU, rhoU, rhoE ],
                 3: [ rho, rhoU, rhoU, rhoU, rhoE ] }[self.dimension]

    def add_default_terms(self):
        self.pde.gamma = self.value('gamma')
        self.pde.R     = self.value('R')
        self.pde.mu    = self.value('mu')
        self.pde.kappa = self.value('kappa')
        self.pde.add_term( name='rhs', type='cf3.sdm.br2.navierstokes_RightHandSide'+str(self.dimension)+'D')
    
    def add_bc(self, name, type, regions, **keyword_args ):
        self.save_bc(name,type,regions,**keyword_args)        
        
        region_comps = [ self.mesh.topology.access_component(str(reg)) for reg in regions ]
          
        if (type == 'Extrapolation' ):
            bc_type = 'cf3.dcm.equations.euler.BCExtrapolation'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'Farfield' ):
            bc_type = 'cf3.dcm.equations.euler.BCFarfield'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            bc.rho = self.value(keyword_args['rho'])
            bc.U = [self.value(v) for v in keyword_args['U'] ]
            bc.p = self.value(keyword_args['p'])
            
        elif (type == 'AdiabaticWall'):
            bc_type = 'cf3.dcm.equations.navierstokes.BCWall'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'SubsonicInlet'):
            bc_type = 'cf3.dcm.equations.euler.BCSubsonicInlet'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            if ('Pt' in keyword_args):
                bc.Pt = self.expression(keyword_args['Pt'])
            if ('Tt' in keyword_args):
                bc.Tt = self.expression(keyword_args['Tt'])
            if ('alpha' in keyword_args):
                bc.alpha = self.expression(keyword_args['alpha'])
            
        elif (type == 'PressureVelocityInlet'):
            bc_type = 'cf3.dcm.equations.euler.BCPressureVelocityInlet'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            if ('p' in keyword_args):
                bc.p = self.expression(keyword_args['p'])
            if ('u' in keyword_args):
                bc.u = self.expression(keyword_args['u'])

        elif (type == 'PressureOutlet'):
            bc_type = 'cf3.dcm.equations.euler.BCPressureOutlet'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            bc.p = self.expression(keyword_args['p'])
            
        else:
            return super(NavierStokes,self).add_bc(name,type,regions)

        return bc

    def print_constants(self):
        log( "Reynolds number =",self.value('Re') )
        log( "density =",self.value('rhoref') )
        log( "Temperature =",self.value('Tref') )
        log( "inlet pressure =",self.value('pref') )
        log( "maximum velocity =",self.value('Uref') )
        log( "sound speed =",self.value('cref') )
        log( "dynamic viscosity =",self.value('mu') )
        log( "kinematic viscosity =",self.value('nu') )
        log( "heat conductivity =",self.value('kappa') )


    def set_post_iteration(self):
        if self.turb_stats_count >= 0:
            self.setup_turbulence_statistics(self.turb_stats_count)
    
    def setup_turbulence_statistics(self,count):
        
        log('Continue turbulence_statistics at count '+str(self.turb_stats_count))
        post_iteration = self.model.get_child('post_iteration')
        if (post_iteration is None):
            post_iteration = self.model.create_component("post_iteration","cf3.common.ActionDirector")

        self.solver.options.post_iteration = post_iteration

        compute_prim_vars = post_iteration.get_child('compute_primitive_variables')
        if (compute_prim_vars is None):
            compute_prim_vars = post_iteration.create_component("compute_primitive_variables",
                "cf3.dcm.equations.navierstokes.ComputePrimitiveVariables")
        compute_prim_vars.gamma = self.value('gamma')
        compute_prim_vars.R = self.value('R')
        compute_prim_vars.solution = self.pde.fields.solution

        self.turb_stats = post_iteration.get_child('compute_turbulence_statistics')
        if (self.turb_stats is None):
            self.turb_stats = post_iteration.create_component('compute_turbulence_statistics',
                'cf3.solver.actions.TurbulenceStatistics')
        self.turb_stats.variable_name = 'U'
        self.turb_stats.region = self.mesh.topology
        self.turb_stats.options.count = count


    def postprocessing(self):

        density = self.pde.fields.get_child('density')
        if( density is None ):
            density = self.pde.fields.create_field(name='density',variables='rho')

        velocity = self.pde.fields.get_child('velocity')
        if( velocity is None ):
            velocity = self.pde.fields.create_field(name='velocity',variables='U[v]')

        pressure = self.pde.fields.get_child('pressure')
        if( pressure is None ):
            pressure = self.pde.fields.create_field(name='pressure',variables='p')

        temperature = self.pde.fields.get_child('temperature')
        if( temperature is None ):
            temperature = self.pde.fields.create_field(name='temperature',variables='T')

        mach = self.pde.fields.get_child('mach')
        if( mach is None ):
            mach = self.pde.fields.create_field(name='mach',variables='M')

        entropy = self.pde.fields.get_child('entropy')
        if( entropy is None ):
            entropy = self.pde.fields.create_field(name='entropy',variables='S')

        total_pressure = self.pde.fields.get_child('total_pressure')
        if( total_pressure is None ):
            total_pressure = self.pde.fields.create_field(name='total_pressure',variables='Pt')

        total_temperature = self.pde.fields.get_child('total_temperature')
        if( total_temperature is None ):
            total_temperature = self.pde.fields.create_field(name='total_temperature',variables='Tt')


        solution = self.pde.fields.solution

        R = self.pde.R
        g = self.pde.gamma

        Sref = self.value('pref/rhoref**gamma')

        for sol,rho,U,p,T,M,S,Pt,Tt in zip(solution,density,velocity,pressure,temperature,mach,entropy,total_pressure,total_temperature):
            rho[0] = sol[0]
            U2 = 0
            for d in range(self.dimension):
                U[d] = sol[1+d]/rho[0]
                U2 += U[d]**2
            p[0] = (g-1)*( sol[3] - 0.5*rho[0]*U2 );
            T[0] = p[0]/(rho[0]*R);
            c = math.sqrt(g*p[0]/rho[0])
            M[0] = math.sqrt(U2)/c
            S[0] = (p[0]/rho[0]**g-Sref)/Sref;
            Pt[0] = p[0]+0.5*rho[0]*U2;
            Tt[0] = T[0]*(1.+((g-1.)/2.))*M[0]**2;

        compute_cfl = self.model.tools.get_child('compute_cfl')
        if ( compute_cfl is None ):
            compute_cfl = self.model.tools.create_component('compute_cfl','cf3.solver.ComputeCFL')

        compute_cfl.wave_speed = self.pde.fields.wave_speed
        compute_cfl.time_step = self.pde.fields.dt
        compute_cfl.execute()

        #self.bdry_layer_regions = ['cylinder']
        bdry_layer_fields = self.pde.bdry_fields
    
        tau_w = bdry_layer_fields.get_child('wall_shear_stress')
        if( tau_w is None ):
            tau_w = bdry_layer_fields.create_field(name='wall_shear_stress',variables='tau_w')
         
        yplus = bdry_layer_fields.get_child('yplus')
        if( yplus is None ):
            yplus = bdry_layer_fields.create_field(name='yplus',variables='yplus')

        ustar = bdry_layer_fields.get_child('ustar')
        if( ustar is None ):
            ustar = bdry_layer_fields.create_field(name='ustar',variables='ustar')
        
        compute_boundary_layer = self.model.tools.get_child('compute_boundary_layer')
        if ( compute_boundary_layer is None ):
            compute_boundary_layer = self.model.tools.create_component('compute_boundary_layer','cf3.dcm.equations.navierstokes.ComputeBoundaryLayer')
        compute_boundary_layer.velocity = self.pde.fields.velocity
        compute_boundary_layer.density = self.pde.fields.density
        compute_boundary_layer.nu = self.value('nu')
        compute_boundary_layer.y0 = self.value('y0')
        compute_boundary_layer.wall_regions = [ self.mesh.topology.access_component(str(reg)) for reg in self.bdry_layer_regions ]
        compute_boundary_layer.yplus = yplus
        compute_boundary_layer.tau = tau_w
        compute_boundary_layer.ustar = ustar
        compute_boundary_layer.execute()
        log("y-plus: between "+str(compute_boundary_layer.yplus_min)+" and "+str(compute_boundary_layer.yplus_max) )

    def __getstate__(self):
        state = super(NavierStokes, self).__getstate__()
        if ( hasattr(self,'turb_stats') ) :
            state['turb_stats_count'] = self.turb_stats.options.count
            log('turb_stats_count = '+str(state['turb_stats_count']))
        state['bdry_layer_regions'] = self.bdry_layer_regions
        return state

    def __setstate__(self,state):
        super(NavierStokes, self).__setstate__(state)
        self.turb_stats_count = state.get('turb_stats_count', -1)
        self.bdry_layer_regions = state.get('bdry_layer_regions',['.'])
