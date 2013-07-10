// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file dcm/erkls_rungekutta/LibERKLS.hpp
/// @author Willem Deconinck
///
/// This file defines the cf3::smd::erkls_rungekutta namespace
/// and the LibERKLS library

#ifndef cf3_dcm_solver_erkls_LibERKLS_hpp
#define cf3_dcm_solver_erkls_LibERKLS_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro dcm_solver_erkls_API
/// @note build system defines COOLFLUID_DCM_EXPLICIT_RUNGEKUTTA_EXPORTS when compiling dcm files
#ifdef COOLFLUID_DCM_SOLVER_ERKLS_EXPORTS
#   define dcm_solver_erkls_API      CF3_EXPORT_API
#   define dcm_solver_erkls_TEMPLATE
#else
#   define dcm_solver_erkls_API      CF3_IMPORT_API
#   define dcm_solver_erkls_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace solver {
/// @brief Explicit time-integration namespace
///
/// @author Willem Deconinck
namespace erkls {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Explicit Runge-Kutta time-integration library
///
/// @author Willem Deconinck
class dcm_solver_erkls_API LibERKLS :
    public cf3::common::Library
{
public:

  /// Constructor
  LibERKLS ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibERKLS() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.dcm.solver.erkls"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "erkls"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements the Explicit Runge-Kutta time-integration method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibERKLS"; }
  
  virtual void initiate();
}; // end LibERKLS

////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_solver_erkls_LibERKLS_hpp
