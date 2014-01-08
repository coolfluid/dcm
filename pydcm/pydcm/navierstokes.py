from pydcm.dcm import DCM
from pydcm.util import log
import math

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

    def __getstate__(self):
        state = super(NavierStokes,self).__getstate__()
        return state
        
    def postprocessing(self):

        density = self.pde.fields.get_child('density')
        if( not density ):
            density = self.pde.fields.create_field(name='density',variables='rho')

        velocity = self.pde.fields.get_child('velocity')
        if( not velocity ):
            velocity = self.pde.fields.create_field(name='velocity',variables='U[v]')

        pressure = self.pde.fields.get_child('pressure')
        if( not pressure ):
            pressure = self.pde.fields.create_field(name='pressure',variables='p')

        temperature = self.pde.fields.get_child('temperature')
        if( not temperature ):
            temperature = self.pde.fields.create_field(name='temperature',variables='T')

        mach = self.pde.fields.get_child('mach')
        if( not mach ):
            mach = self.pde.fields.create_field(name='mach',variables='M')

        entropy = self.pde.fields.get_child('entropy')
        if( not entropy ):
            entropy = self.pde.fields.create_field(name='entropy',variables='S')

        total_pressure = self.pde.fields.get_child('total_pressure')
        if( not total_pressure ):
            total_pressure = self.pde.fields.create_field(name='total_pressure',variables='Pt')

        total_temperature = self.pde.fields.get_child('total_temperature')
        if( not total_temperature ):
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
