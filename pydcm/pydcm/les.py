from pydcm.navierstokes import NavierStokes

class LES(NavierStokes):
    
    def __init__(self,name='les',dimension=0):
        super(LES, self).__init__(name,dimension)
        self.pde_type = 'cf3.dcm.equations.les.LES'
           
    def add_default_terms(self):
        self.pde.gamma = self.value('gamma')
        self.pde.R     = self.value('R')
        self.pde.mu    = self.value('mu')
        self.pde.kappa = self.value('kappa')
        self.pde.sfs_model = 'WALE'
        self.pde.riemann_solver = 'Rusanov'
        self.pde.add_term( name='rhs', type='cf3.sdm.br2.les_RightHandSide'+str(self.dimension)+'D')

    def postprocessing(self):
        super(LES, self).postprocessing()
        
        sfs_kinetic_energy = self.pde.fields.get_child('sfs_kinetic_energy')
        if( sfs_kinetic_energy is None ):
            sfs_kinetic_energy = self.pde.fields.create_field(name='sfs_kinetic_energy',variables='k_sfs')

        sfs_heat_conduction = self.pde.fields.get_child('sfs_heat_conduction')
        if( sfs_heat_conduction is None ):
            sfs_heat_conduction = self.pde.fields.create_field(name='sfs_heat_conduction',variables='kappaT')

        sfs_viscosity = self.pde.fields.get_child('sfs_viscosity')
        if( sfs_viscosity is None ):
            sfs_viscosity = self.pde.fields.create_field(name='sfs_viscosity',variables='nuT')
        
        compute_subfilterscale = self.model.tools.get_child('compute_subfilterscale')
        if ( compute_subfilterscale is None ):
            compute_subfilterscale = self.model.tools.create_component('compute_subfilterscale','cf3.dcm.equations.les.ComputeSubFilterScale')
        # wale = compute_subfilterscale.create_component('wale','cf3.dcm.equations.les.WALE2D')
        compute_subfilterscale.velocity = self.pde.fields.velocity
        compute_subfilterscale.density = self.pde.fields.density
        compute_subfilterscale.sfs_model = self.pde.access_component('rhs_computer/rhs/term/sfs_model')
        # compute_subfilterscale.sfs_model = wale
        compute_subfilterscale.sfs_kinetic_energy = sfs_kinetic_energy
        compute_subfilterscale.sfs_heat_conduction = sfs_heat_conduction
        compute_subfilterscale.sfs_viscosity = sfs_viscosity
        compute_subfilterscale.execute()
        