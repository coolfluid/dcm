// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file dcm/erk_rungekutta/LibERK.hpp
/// @author Willem Deconinck
///
/// This file defines the cf3::smd::erk_rungekutta namespace
/// and the LibERK library

#ifndef cf3_dcm_solver_erk_LibERK_hpp
#define cf3_dcm_solver_erk_LibERK_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro dcm_solver_erk_API
/// @note build system defines COOLFLUID_DCM_SOLVER_ERK_EXPORTS when compiling dcm files
#ifdef COOLFLUID_DCM_SOLVER_ERK_EXPORTS
#   define dcm_solver_erk_API      CF3_EXPORT_API
#   define dcm_solver_erk_TEMPLATE
#else
#   define dcm_solver_erk_API      CF3_IMPORT_API
#   define dcm_solver_erk_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace solver {
/// @brief ERK time-integration namespace
///
/// @author Willem Deconinck
namespace erk {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the ERK Runge-Kutta time-integration library
///
/// @author Willem Deconinck
class dcm_solver_erk_API LibERK :
    public cf3::common::Library
{
public:

  /// Constructor
  LibERK ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibERK() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.dcm.solver.erk"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "erk"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements the ERK Runge-Kutta time-integration method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibERK"; }
  
  virtual void initiate();
}; // end LibERK

////////////////////////////////////////////////////////////////////////////////

} // erk
} // solver
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_solver_erk_LibERK_hpp
