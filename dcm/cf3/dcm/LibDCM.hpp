// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_LibDCM_hpp
#define cf3_dcm_LibDCM_hpp

#define SANDBOX
//#undef SANDBOX
////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro dcm_API
/// @note build system defines COOLFLUID_DCM_EXPORTS when compiling dcm files
#ifdef COOLFLUID_DCM_EXPORTS
#   define dcm_API      CF3_EXPORT_API
#   define dcm_TEMPLATE
#else
#   define dcm_API      CF3_IMPORT_API
#   define dcm_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {

/// @brief Spectral Finite Difference Method namespace
///
/// The Spectral Finite Difference Method is a high-order method
/// for solving systems of partial differential equations.
/// @author Willem Deconinck
namespace dcm {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Spectral Finite Difference Core library
///
/// This library implements Core components to construct a Spectral Finite Difference Solver.";
/// @author Willem Deconinck
class dcm_API LibDCM :
    public cf3::common::Library
{
public:

  
  

  /// Constructor
  LibDCM ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibDCM() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.dcm"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "dcm"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements Core components to construct a Spectral Finite Difference Solver.";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibDCM"; }
  
  virtual void initiate();
}; // end LibDCM

////////////////////////////////////////////////////////////////////////////////

} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_LibDCM_hpp
