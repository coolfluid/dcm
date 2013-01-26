// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_advectiondiffusion_LibAdvectionDiffusion_hpp
#define cf3_sdm_equations_advectiondiffusion_LibAdvectionDiffusion_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Library.hpp"

////////////////////////////////////////////////////////////////////////////////

/// Define the macro advectiondiffusion_API
/// @note build system defines COOLFLUID_SDM_EQUATIONS_ADVECTIONDIFFUSION_EXPORTS when compiling advectiondiffusion files
#ifdef COOLFLUID_SDM_EQUATIONS_ADVECTIONDIFFUSION_EXPORTS
#   define sdm_equations_advectiondiffusion_API      CF3_EXPORT_API
#   define sdm_equations_advectiondiffusion_TEMPLATE
#else
#   define sdm_equations_advectiondiffusion_API      CF3_IMPORT_API
#   define sdm_equations_advectiondiffusion_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

/// Class defines the advectiondiffusion library
class sdm_equations_advectiondiffusion_API LibAdvectionDiffusion : public common::Library
{
public:

  /// Constructor
  LibAdvectionDiffusion ( const std::string& name) : common::Library(name) { }

  virtual ~LibAdvectionDiffusion() { }

public: // functions

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.sdm.equations.advectiondiffusion"; }

  /// Static function that returns the library name.
  /// Must be implemented for Library registration
  /// @return name of the library
  static std::string library_name() { return "advectiondiffusion"; }

  /// Static function that returns the description of the library.
  /// Must be implemented for Library registration
  /// @return description of the library

  static std::string library_description()
  {
    return "This library implements the Advection-Diffusion equation for the SD method";
  }

  /// Gets the Class name
  static std::string type_name() { return "LibAdvectionDiffusion"; }

}; // end LibAdvectionDiffusion

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_advectiondiffusion_LibAdvectionDiffusion_hpp

