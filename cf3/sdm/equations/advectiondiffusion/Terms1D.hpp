// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_advectiondiffusion_Terms1D_hpp
#define cf3_sdm_equations_advectiondiffusion_Terms1D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/math/Checks.hpp"
#include "cf3/solver/Term.hpp"
#include "cf3/sdm/equations/advectiondiffusion/LibAdvectionDiffusion.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

  class sdm_equations_advectiondiffusion_API Terms1D : public solver::TermBase< 1 /*dim*/, 1 /*eqs*/, 1 /*vars*/, 1/*grads*/ >
{
public: 

  /// @brief Constructor
  Terms1D( const std::string& name );
  
  /// @brief Destructor
  virtual ~Terms1D() {}
  
  static std::string type_name() { return "Terms1D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_DIFFUSION  = true };

  struct DATA {
    Real q;       // solution variable
    Real grad_q;  // solution gradient
    Real a;       // advection speed
    Real mu;      // diffusion coefficient
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
                                       RowVector_NEQS& flux, Real& wave_speed );

  static void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                                    RowVector_NEQS& flux, Real& wave_speed );

  static void compute_diffusive_flux( const DATA& p, const ColVector_NDIM& normal,
                                      RowVector_NEQS& flux, Real& wave_speed );

private: // configuration

  Real m_a;
  Real m_mu;
};

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_advectiondiffusion_Terms1D_hpp
