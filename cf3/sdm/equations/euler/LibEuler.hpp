// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_euler_LibEuler_hpp
#define cf3_sdm_equations_euler_LibEuler_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro euler_API
/// @note build system defines COOLFLUID_SDM_EULER_EXPORTS when compiling euler files
#ifdef COOLFLUID_SDM_EQUATIONS_EULER_EXPORTS
#   define sdm_equations_euler_API      CF3_EXPORT_API
#   define sdm_equations_euler_TEMPLATE
#else
#   define sdm_equations_euler_API      CF3_IMPORT_API
#   define sdm_equations_euler_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

/// Class defines the euler library
class sdm_equations_euler_API LibEuler : public common::Library
{
public:

  /// Constructor
  LibEuler ( const std::string& name) : common::Library(name) { }

  virtual ~LibEuler() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.sdm.equations.euler"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "euler"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements the Euler equations for the SD method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibEuler"; }

}; // end LibEuler

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_euler_LibEuler_hpp

