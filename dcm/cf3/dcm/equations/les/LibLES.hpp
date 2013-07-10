// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_les_LibLES_hpp
#define cf3_dcm_equations_les_LibLES_hpp

////////////////////////////////////////////////////////////////////////////////

#include "common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro les_API
/// @note build system defines COOLFLUID_dcm_LES_EXPORTS when compiling les files
#ifdef COOLFLUID_DCM_EQUATIONS_LES_EXPORTS
#   define dcm_equations_les_API      CF3_EXPORT_API
#   define dcm_equations_les_TEMPLATE
#else
#   define dcm_equations_les_API      CF3_IMPORT_API
#   define dcm_equations_les_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

/// Class defines the les library
class dcm_equations_les_API LibLES : public common::Library
{
public:

  /// Constructor
  LibLES ( const std::string& name) : common::Library(name) { }

  virtual ~LibLES() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.dcm.equations.les"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "les"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements les equations";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibLES"; }

}; // end LibLES

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_les_LibLES_hpp

