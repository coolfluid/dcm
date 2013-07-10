// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_les_Smagorinsky2D_hpp
#define cf3_dcm_equations_les_Smagorinsky2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/Term.hpp"
#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/physics/navierstokes/navierstokes2d/Functions.hpp"
#include "cf3/dcm/equations/navierstokes/RightHandSide2D.hpp"
#include "cf3/dcm/equations/les/LibLES.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
  namespace mesh {
    class Field;
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_les_API Smagorinsky2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 5 /*vars*/, 3/*grads*/ >
{
public: 

  /// @brief Constructor
  Smagorinsky2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~Smagorinsky2D() {}
  
  static std::string type_name() { return "LESSmagorinsky2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_DIFFUSION  = true };

  struct DATA: physics::navierstokes::navierstokes2d::Data
  {
    
    Real PrT;               ///< Turbulent Prandtl number      ( 0.7 - 0.9 )
    Real Cs;                ///< Smagorinsky constant          ( 0.1 < Cs < 0.24 )
    Real Cv;                ///< Deardorff/Yoshizawa constant  ( 0.094 )
    Real k_sfs;
    Real nuT;
    Real muT;
    Real kappaT;
    Real Delta;             ///< Turbulent mixing length scale
  };



public: // Variable and PhysData computation
    
  /// @brief Compute variables and gradients in a given element point
  ///
  /// The interpolation and gradient reconstructions, as well as
  void get_variables( const mesh::Space& space,
                      const Uint elem_idx,
                      const ColVector_NDIM& coords,
                      const mesh::ReconstructPoint& interpolation,
                      const std::vector<mesh::ReconstructPoint>& gradient,
                      const Matrix_NDIMxNDIM& jacobian,
                      const Matrix_NDIMxNDIM& jacobian_inverse,
                      const Real& jacobian_determinant,
                      RowVector_NVAR& vars,
                      RowVector_NGRAD& gradvars,
                      Matrix_NDIMxNGRAD& gradvars_grad );
  
  void get_bdry_variables( const mesh::Space& space,
                           const Uint elem_idx,
                           const ColVector_NDIM& coords,
                           const mesh::ReconstructPoint& interpolation,
                           const std::vector<mesh::ReconstructPoint>& gradient,
                           const Matrix_NDIMxNDIM& jacobian,
                           const Matrix_NDIMxNDIM& jacobian_inverse,
                           const Real& jacobian_determinant,
                           RowVector_NVAR& vars,
                           RowVector_NGRAD& gradvars,
                           Matrix_NDIMxNGRAD& gradvars_grad );

  /// @brief Set constants in the data
  void set_phys_data_constants( DATA& phys_data );

  /// @brief Compute the data from computed variables and gradients
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Flux computations

  static void compute_convective_flux( const DATA& p, const ColVector_NDIM& normal,
                                       RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_convective_flux(p,normal,flux,wave_speed);
  }

  void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                             RowVector_NEQS& flux, Real& wave_speed )
  {
    m_riemann_solver->compute_riemann_flux(left,right,normal,flux,wave_speed);
  }

  static void compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed);


private:

  void config_riemann_solver();

private: // configuration

  Handle< solver::RiemannSolver<cf3::physics::euler::euler2d::Data,NDIM,NEQS> > m_riemann_solver;
  std::string m_riemann_solver_type;
  Real m_gamma;
  Real m_R;
  Real m_kappa;
  Real m_mu;
  Real m_PrT;
  Real m_Cs;
  Real m_Cv;


  Real _rho;
  Real _u;
  Real _v;
  Real _E;
  Real _U2;
  Real _p;
  Real _T;
  ColVector_NDIM _grad_rho;
  ColVector_NDIM _grad_rhou;
  ColVector_NDIM _grad_rhov;
  ColVector_NDIM _grad_rhoE;
  ColVector_NDIM _grad_u;
  ColVector_NDIM _grad_v;
  ColVector_NDIM _grad_T;

};

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_les_Smagorinsky2D_hpp
