// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_core_LibCore_hpp
#define cf3_sdm_core_LibCore_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro sdm_core_API
/// @note build system defines COOLFLUID_SDM_EXPORTS when compiling sdm files
#ifdef COOLFLUID_SDM_CORE_EXPORTS
#   define sdm_core_API      CF3_EXPORT_API
#   define sdm_core_TEMPLATE
#else
#   define sdm_core_API      CF3_IMPORT_API
#   define sdm_core_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {

/// @brief Spectral Finite Difference Method namespace
///
/// The Spectral Finite Difference Method is a high-order method
/// for solving systems of partial differential equations.
/// @author Willem Deconinck
namespace sdm {
  
/// @brief Core Spectral Finite Difference Method namespace
///
/// The Spectral Finite Difference Method is a high-order method
/// for solving systems of partial differential equations.
/// @author Willem Deconinck
namespace core {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Spectral Finite Difference Core library
///
/// This library implements Core components to construct a Spectral Finite Difference Solver.";
/// @author Willem Deconinck
class sdm_core_API LibCore :
    public cf3::common::Library
{
public:

  /// Constructor
  LibCore ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibCore() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.sdm.core"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "core"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements Core components to construct a Spectral Finite Difference Solver.";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibCore"; }
  
  virtual void initiate();
  
}; // end LibCore

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_core_LibCore_hpp
