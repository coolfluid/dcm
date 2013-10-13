from pydcm.dcm import DCM

class NavierStokes(DCM):
    
    def __init__(self,name='navierstokes',dimension=0):
        super(NavierStokes, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.navierstokes.NavierStokes'
                        
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

    def __getstate__(self):
        state = super(NavierStokes,self).__getstate__()
        return state
        

class NavierStokesParam(object):
    
    def __init__(self):
        # gas constants for AIR
        self.gamma = 1.4
        self.R = 287.05
        self.Pr = 0.72
        self.kappa = 2.601e-2
        self.mu = 1.806e-5
        self.rho = 1.2041
        self.p = 1.e5
        self.M = 0.
        self.u = 0.
        self.v = 0.
        self.w = 0.


    def compute_conservative_state(self,state):
        # flow parameters
        if len(state) == 3:
            self.rho = state[0]
            self.u = state[1]/rho
            self.E = state[2]/rho
        elif len(state) == 4:
            self.rho = state[0]
            self.u = state[1]/rho
            self.v = state[2]/rho
            self.E = state[3]/rho
        elif len(state) == 5:
            self.rho = state[0]
            self.u = state[1]/rho
            self.v = state[2]/rho
            self.w = state[3]/rho
            self.E = state[4]/rho
        self.U  = sqrt(self.u**2+self.v**2+self.w**2)
        self.p  = (self.gamma-1)*(self.rho*self.E - 0.5*self.rho*self.U**2);
        self.c  = sqrt(self.gamma*self.p/self.rho)
        self.T  = self.p/(self.rho*self.R);
        self.M  = self.U/self.c;
        self.Pt = self.p+0.5*self.rho*self.U**2;
        self.Tt = self.T*(1.+((self.gamma-1)/2))*self.M**2;
        self.S  = self.p/(abs(self.rho)**self.gamma);
        self.T  = self.p / (self.rho * self.R)
        self.nu = self.mu/self.rho
        
        
        
        
