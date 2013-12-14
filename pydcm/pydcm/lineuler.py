from pydcm.dcm import DCM
from math import sqrt

class LinearizedEuler(DCM):
    
    def __init__(self,name='lineuler',dimension=0):
        super(LinearizedEuler, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.lineuler.LinEuler'
                        
    def ref_solution(self):
        rho0 = self.value('rho0')
        rho  = rho0*0.01
        p = self.value('p0')*0.01
        c0 = sqrt(self.value('gamma')*self.value('p0')/self.value('rho0'))
        U = 0.01*c0
        rhoU = rho*U
        return { 1: [ rho, rho0*U, p ],
                 2: [ rho, rho0*U, rho0*U, p ],
                 3: [ rho, rho0*U, rho0*U, rho0*U, p ] }[self.dimension]

    def add_default_terms(self):
        self.pde.gamma = self.value('gamma')
        self.pde.add_term( name='rhs', type='cf3.sdm.br2.lineuler_RightHandSide'+str(self.dimension)+'D')
    
    def add_bc(self, name, type, regions, **keyword_args ):
        self.save_bc(name,type,regions,**keyword_args)        
        
        region_comps = [ self.mesh.topology.access_component(str(reg)) for reg in regions ]
          
        if (type == 'Extrapolation' ):
            bc_type = 'cf3.dcm.equations.lineuler.BCExtrapolation'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'Farfield' ):
            bc_type = 'cf3.dcm.equations.lineuler.BCFarfield'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            
        elif (type == 'Mirror' or type == 'Reflection'):
            bc_type = 'cf3.dcm.equations.lineuler.BCMirror'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'Thompson' ):
            bc_type = 'cf3.sdmx.equations.lineuler.BCThompsonUniform'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        elif (type == 'NonReflecting' ):
            bc_type = 'cf3.sdmx.equations.lineuler.BCNonReflecting'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )

        else:
            return super(LinearizedEuler,self).add_bc(name,type,regions)

        return bc
        
    def init_background(self, *functions):
        self.pde.fields.print_tree()
        
        for bg_field in [ self.pde.fields.background, self.pde.bdry_fields.bdry_background ]:
            self.model.tools.init_field.init_field(
              field=bg_field,
              functions=[self.expression(func) for func in functions] )
                       
        gradient_computer = self.model.tools.create_component('gradient_computer','cf3.dcm.tools.ComputeFieldGradientBR2')
        gradient_computer.field = self.pde.fields.background
        gradient_computer.field_gradient = self.pde.fields.background_gradient
        gradient_computer.execute()
