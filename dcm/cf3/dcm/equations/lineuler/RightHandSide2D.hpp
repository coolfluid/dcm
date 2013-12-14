// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_RightHandSide2D_hpp
#define cf3_dcm_equations_lineuler_RightHandSide2D_hpp

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

class dcm_equations_lineuler_API RightHandSide2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 16 /*vars*/, 0/*grads*/ >
{
public: 
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// @brief Constructor
  RightHandSide2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~RightHandSide2D() {}
  
  static std::string type_name() { return "RightHandSide2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_SOURCE     = true };

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

  /// @brief Source term, coming from inhomogeneous mean flow
  ///
  /// The source term H is defined as
  ///
  ///      [                                 0                                  ]
  ///  H = [     (rho' u0 + rho0 u') du0/dx  +  (rho' v0 + rho0 v') du0/dy      ]
  ///      [     (rho' u0 + rho0 u') dv0/dx  +  (rho' v0 + rho0 v') dv0/dy      ]
  ///      [ (gamma-1) p' (du0/dx + dv0/dy) - (gamma-1) (u' dp0/dx + v' dp0/dy) ]
  ///
  /// with u0 = U0[XX]
  /// with v0 = U0[YY]
  /// with u' = U[XX] - U0[XX]
  /// with v' = U[YY] - U0[YY]
  void compute_source( const DATA& p, RowVector_NEQS& source );

private: // configuration

  Real m_gamma;
  Handle<mesh::Field> m_background;
  Handle<mesh::Field> m_background_gradient;
  Handle<mesh::Field> m_bdry_background;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_RightHandSide2D_hpp
