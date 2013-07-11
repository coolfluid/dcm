// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_br2_LibBR2_hpp
#define cf3_sdm_br2_LibBR2_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro sdm_br2_API
/// @note build system defines COOLFLUID_SDM_BR2_EXPORTS when compiling
#ifdef COOLFLUID_SDM_BR2_EXPORTS
#   define sdm_br2_API      CF3_EXPORT_API
#   define sdm_br2_TEMPLATE
#else
#   define sdm_br2_API      CF3_IMPORT_API
#   define sdm_br2_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {

/// @brief Spectral Finite Difference Method namespace
///
/// The Spectral Finite Difference Method is a high-order method
/// for solving systems of partial differential equations.
/// @author Willem Deconinck
namespace sdm {

/// @brief BR2 namespace
///
/// Library containing term computers of the Spectral Difference Method.
/// Diffusive terms are computed using the second approach of Bassi-Rebay (BR2)
/// @author Willem Deconinck
namespace br2 {

////////////////////////////////////////////////////////////////////////////////

/// @brief Defines the Spectral Finite Difference library with BR2 scheme
/// @author Willem Deconinck
class sdm_br2_API LibBR2 :
    public cf3::common::Library
{
public:

  /// Constructor
  LibBR2 ( const std::string& name) : cf3::common::Library(name) { }

  virtual ~LibBR2() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.sdm.br2"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "br2"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements SDM with BR2 scheme";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibBR2"; }
  
  virtual void initiate();
}; // end LibBR2

////////////////////////////////////////////////////////////////////////////////

} // br2
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_br2_LibBR2_hpp
