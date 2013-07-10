// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_advectiondiffusion_AdvectionDiffusion1D_hpp
#define cf3_dcm_equations_advectiondiffusion_AdvectionDiffusion1D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/PDE.hpp"
#include "cf3/dcm/equations/advectiondiffusion/LibAdvectionDiffusion.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

/// @brief advectiondiffusion 1D physics
///
/// This component assembles the terms, term-computers,
/// and configuration options
/// to solve the right-hand-side of the advectiondiffusion-equations in
/// the form  dQ/dt = RHS
/// @author Willem Deconinck
class dcm_equations_advectiondiffusion_API AdvectionDiffusion1D : public solver::PDE {

public: // functions

  /// Contructor
  /// @param name of the component
  AdvectionDiffusion1D ( const std::string& name );

  /// Virtual destructor
  virtual ~AdvectionDiffusion1D();

  /// Get the class name
  static std::string type_name () { return "AdvectionDiffusion1D"; }

};

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_advectiondiffusion_AdvectionDiffusion1D_hpp
