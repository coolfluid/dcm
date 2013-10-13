from pydcm.dcm import DCM

class AdvectionDiffusion(DCM):
    
    def __init__(self,name='advection_diffusion',dimension=0):
        super(AdvectionDiffusion, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.advectiondiffusion.AdvectionDiffusion'
                        
    def ref_solution(self):
        if not self.defined('Q_ref'): self.define(Q_ref=1.)            
        return [self.value('Q_ref')]

    def add_default_terms(self):
        self.pde.a  = self.value('a')
        self.pde.mu = self.value('mu')
        self.pde.add_term( name='rhs', type='cf3.sdm.br2.advectiondiffusion_RightHandSide'+str(self.dimension)+'D')        
        
    def add_bc(self, name, type, regions, **keyword_args ):
        self.save_bc(name,type,regions,**keyword_args)        
        region_comps = [ self.mesh.topology.access_component(str(reg)) for reg in regions ]
    
        if (type == 'Extrapolation' ):
            bc_type = 'cf3.dcm.equations.advectiondiffusion.BCExtrapolation'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
        elif (type == 'Dirichlet' ):
            bc_type = 'cf3.dcm.equations.advectiondiffusion.BCDirichlet'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            if 'Q' in keyword_args: bc.Q = keyword_args['Q']
        elif (type == 'Neumann' ):
            bc_type = 'cf3.dcm.equations.advectiondiffusion.BCNeumann'+str(self.dimension)+'D'
            bc = self.pde.add_bc( name=name, type=bc_type, regions=region_comps )
            if 'dQdn' in keyword_args: bc.dQdn = keyword_args['dQdn']            
        else:
            return super(AdvectionDiffusion,self).add_bc(name,type,regions)
            
    