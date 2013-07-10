// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file plugins/sdm/cf3/sdm/solver/optim_erkls/Types.hpp
/// @author Willem Deconinck
/// @author Matteo Parsani
///
/// This file defines the coefficients of some low-storage explicit Runge-Kutta
/// schemes optimized for the spectral difference method
/// A integration scheme of order X is optimised for a SDM of the same order ( P=order-1 )
///
/// Derivation and coefficients are taken from:
///   Parsani, M., Ketcheson, D. I., & Deconinck, W. (2013). 
///   Optimized Explicit Runge-Kutta Schemes for the Spectral Difference Method Applied to Wave Propagation Problems. 
///   SIAM J. Scientific Computing (SIAMSC) 35(2). doi:10.1137/120885899

#ifndef cf3_sdm_solver_optim_erkls_Types_hpp
#define cf3_sdm_solver_optim_erkls_Types_hpp

#include "cf3/sdm/solver/optim_erkls/LibOptimERKLS.hpp"
#include "cf3/dcm/solver/erkls/ThreeSstar.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace solver {
namespace optim_erkls {

/////////////////////////////////////////////////////////////////////////////////////

// Macro to declare ThreeSstarCoeffs::Coefficients
// in the namespace "erk_coeffs", with name
// ERK_stage_order
// Also a ThreeSstarT<ERK_stage_order> is defined in current namespace
// with name ERK_stage_order as well.
#define DECLARE_ERK(stage,order) \
namespace erk_coeffs { \
  struct ERK_##stage##_##order : dcm::solver::erkls::ThreeSstarCoeffs::Coefficients \
  { \
    static std::string name() { return "ERK_"+std::string(#stage)+"_"+std::string(#order); } \
    ERK_##stage##_##order(); \
  }; \
} \
typedef dcm::solver::erkls::ThreeSstarT< erk_coeffs::ERK_##stage##_##order  >  ERK_##stage##_##order;

/////////////////////////////////////////////////////////////////////////////////////

// Second order schemes, optimized for Second order Spectral Difference (P1)
DECLARE_ERK(3,2);   // stable CFL = 0.587624777108
DECLARE_ERK(6,2);   // stable CFL = 1.251104269176
DECLARE_ERK(8,2);   // stable CFL = 1.677445545792

// Third order schemes, optimized for Third order Spectral Difference (P2)
DECLARE_ERK(5,3);   // stable CFL = 0.453596841544
DECLARE_ERK(8,3);   // stable CFL = 0.808357894420
DECLARE_ERK(17,3);  // stable CFL = 1.82207767386

// Fourth order schemes, optimized for Fourth order Spectral Difference (P3)
DECLARE_ERK(9,4);   // stable CFL = 0.512795113027
DECLARE_ERK(14,4);  // stable CFL = 0.882074842229
DECLARE_ERK(18,4);  // stable CFL = 1.17418695241

// Fifth order schemes, optimized for Fifth order Spectral Difference (P4)
DECLARE_ERK(10,5);  // stable CFL = 0.361642148346
DECLARE_ERK(20,5);  // stable CFL = 0.843896716833

/////////////////////////////////////////////////////////////////////////////////////

// Undefine the macro again as it no longer needed
#undef DECLARE_ERK

} // optim_erkls
} // solver
} // sdm
} // cf3

#endif // cf3_sdm_solver_optim_erkls_Types_hpp
