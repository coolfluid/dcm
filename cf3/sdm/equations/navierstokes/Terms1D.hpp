// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_navierstokes_Terms1D_hpp
#define cf3_sdm_equations_navierstokes_Terms1D_hpp

////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "cf3/solver/Term.hpp"
#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/physics/euler/euler1d/Functions.hpp"
#include "cf3/physics/navierstokes/navierstokes1d/Functions.hpp"
#include "cf3/sdm/equations/navierstokes/LibNavierStokes.hpp"
namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_navierstokes_API Terms1D : public solver::TermBase< 1 /*dim*/, 3 /*eqs*/, 3 /*vars*/, 2/*grads*/ >
{
public: 

  /// @brief Constructor
  Terms1D( const std::string& name );
  
  /// @brief Destructor
  virtual ~Terms1D() {}
  
  static std::string type_name() { return "NSTerms1D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_DIFFUSION = true };

  typedef physics::navierstokes::navierstokes1d::Data DATA;

public: // Variable and PhysData computation
    
  /// @brief Compute variables and gradients in a given element point
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
    physics::euler::euler1d::compute_convective_flux(p,normal,flux,wave_speed);
    std::cout << "Fc = " << flux << std::endl;
  }

  void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                             RowVector_NEQS& flux, Real& wave_speed )
  {
    m_riemann_solver->compute_riemann_flux(left,right,normal,flux,wave_speed);
    std::cout << "Fc = " << flux << std::endl;
  }

  static void compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed);


private:

  void config_riemann_solver();

private: // configuration

  Handle< solver::RiemannSolver<Terms1D> > m_riemann_solver;
  std::string m_riemann_solver_type;
  Real m_gamma;
  Real m_R;
  Real m_kappa;
  Real m_mu;
};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_navierstokes_Terms1D_hpp
