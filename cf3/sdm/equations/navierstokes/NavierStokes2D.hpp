// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_navierstokes_NavierStokes2D_hpp
#define cf3_sdm_equations_navierstokes_NavierStokes2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/PDE.hpp"
#include "cf3/sdm/equations/navierstokes/LibNavierStokes.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

/// @brief NavierStokes 2D physics
///
/// This component assembles the terms, term-computers,
/// and configuration options
/// to solve the right-hand-side of the NavierStokes-equations in
/// the form  dQ/dt = RHS
/// @author Willem Deconinck
class sdm_equations_navierstokes_API NavierStokes2D : public solver::PDE {

public: // functions

  /// Contructor
  /// @param name of the component
  NavierStokes2D ( const std::string& name );

  /// Virtual destructor
  virtual ~NavierStokes2D();

  /// Get the class name
  static std::string type_name () { return "NavierStokes2D"; }

  virtual std::string solution_variables() const;

};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_navierstokes_NavierStokes2D_hpp
