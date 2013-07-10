// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file dcm/optim_erkls_rungekutta/LibOptimERKLS.hpp
/// @author Willem Deconinck
///
/// This file defines the cf3::smd::optim_erkls_rungekutta namespace
/// and the LibOptimERKLS library

#ifndef cf3_sdm_solver_optim_erkls_LibOptimERKLS_hpp
#define cf3_sdm_solver_optim_erkls_LibOptimERKLS_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro sdm_solver_optim_erkls_API
/// @note build system defines COOLFLUID_SDM_EXPLICIT_RUNGEKUTTA_EXPORTS when compiling sdm files
#ifdef COOLFLUID_SDMX_SOLVER_OPTIM_ERKLS_EXPORTS
#   define sdm_solver_optim_erkls_API      CF3_EXPORT_API
#   define sdm_solver_optim_erkls_TEMPLATE
#else
#   define sdm_solver_optim_erkls_API      CF3_IMPORT_API
#   define sdm_solver_optim_erkls_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace solver {
/// @brief Explicit time-integration namespace
///
/// @author Willem Deconinck
namespace optim_erkls {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Explicit Runge-Kutta time-integration library
///
/// @author Willem Deconinck
class sdm_solver_optim_erkls_API LibOptimERKLS :
    public cf3::common::Library
{
public:

  /// Constructor
  LibOptimERKLS ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibOptimERKLS() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.sdm.solver.optim_erkls"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "optim_erkls"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements the Explicit Runge-Kutta time-integration method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibOptimERKLS"; }
  
  virtual void initiate();
}; // end LibOptimERKLS

////////////////////////////////////////////////////////////////////////////////

} // optim_erkls
} // solver
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_solver_optim_erkls_LibOptimERKLS_hpp
