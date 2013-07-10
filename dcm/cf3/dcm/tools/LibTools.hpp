// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_tools_LibTools_hpp
#define cf3_dcm_tools_LibTools_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro dcm_API
/// @note build system defines COOLFLUID_DCM_EXPORTS when compiling dcm files
#ifdef COOLFLUID_DCM_TOOLS_EXPORTS
#   define dcm_tools_API      CF3_EXPORT_API
#   define dcm_tools_TEMPLATE
#else
#   define dcm_tools_API      CF3_IMPORT_API
#   define dcm_tools_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {

/// @brief Spectral Finite Difference Tools namespace
///
/// The Spectral Finite Difference Method is a high-order method
/// for solving systems of partial differential equations.
/// This namespace contains useful tools
/// @author Willem Deconinck
namespace tools {
  
////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Spectral Finite Difference Core library
///
/// This library implements useful tool components for the Spectral Difference method";
/// @author Willem Deconinck
class dcm_tools_API LibTools :
    public cf3::common::Library
{
public:

  /// Constructor
  LibTools ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibTools() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.dcm.tools"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "tools"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements useful components specific to the Spectral Difference method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibTools"; }
  
  virtual void initiate();
}; // end LibTools

////////////////////////////////////////////////////////////////////////////////

} // tools
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_tools_LibTools_hpp
