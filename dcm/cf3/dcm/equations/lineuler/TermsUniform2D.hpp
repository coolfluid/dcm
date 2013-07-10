// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_TermsUniform2D_hpp
#define cf3_dcm_equations_lineuler_TermsUniform2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/Term.hpp"
#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/physics/lineuler/lineuler2d/Functions.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API TermsUniform2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 4 /*vars*/, 0/*grads*/ >
{
public: 

  /// @brief Constructor
  TermsUniform2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~TermsUniform2D() {}
  
  static std::string type_name() { return "TermsUniform2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };

  typedef physics::lineuler::lineuler2d::Data DATA;

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
    physics::lineuler::lineuler2d::compute_convective_flux(p,normal,flux,wave_speed);
  }

  void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                             RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::lineuler::lineuler2d::compute_cir_flux(left,right,normal,flux,wave_speed);
  }

private: // configuration

  Real m_gamma;
  std::vector<Real> m_U0;
  Real m_rho0;
  Real m_p0;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_TermsUniform2D_hpp
