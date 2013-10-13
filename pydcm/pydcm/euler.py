from pydcm.dcm import DCM

class Euler(DCM):
    
    def __init__(self,name='euler',dimension=0):
        super(Euler, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.euler.Euler'
                        
    def ref_solution(self):
        rho  = self.value('rho_ref')
        U = self.value('U_ref')
        p = self.value('p_ref')
        rhoU = rho*U
        rhoE = p/(self.value('gamma')-1.)+0.5*rho*U**2
        return { 1: [ rho, rhoU, rhoE ],
                 2: [ rho, rhoU, rhoU, rhoE ],
                 3: [ rho, rhoU, rhoU, rhoU, rhoE ] }[self.dimension]

    def add_default_terms(self):
        self.pde.gamma = self.value('gamma')
        self.pde.R     = self.value('R')
        self.pde.add_term( name='rhs', type='cf3.sdm.br2.euler_RightHandSide'+str(self.dimension)+'D')
    
    def add_bc(self, name, type, regions, **keyword_args ):
        self.save_bc(name,type,regions,**keyword_args)        
        
        region_comps = [ self.mesh.topology.access_component(str(reg)) for reg in regions ]
          
        if (type == 'Extrapolation' ):
            bc_type = 'cf3.dcm.equations.euler.BCExtrapolation'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'Farfield' ):
            bc_type = 'cf3.dcm.equations.euler.BCFarfield'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            
        elif (type == 'Mirror' or type == 'Reflection'):
            bc_type = 'cf3.dcm.equations.euler.BCMirror'+str(self.dimension)+'D'
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
                bc.U = self.expression(keyword_args['u'])

        elif (type == 'PressureOutlet'):
            bc_type = 'cf3.dcm.equations.euler.BCPressureOutlet'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            bc.p = self.expression(keyword_args['p'])
            
        else:
            return super(Euler,self).add_bc(name,type,regions)

        return bc

    def print_constants(self):
        log( "density =",self.value('rhoref') )
        log( "Temperature =",self.value('Tref') )
        log( "inlet pressure =",self.value('pref') )
        log( "maximum velocity =",self.value('Uref') )
        log( "sound speed =",self.value('cref') )
        
        