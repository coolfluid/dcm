// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_lineuler_LinEulerUniform2D_hpp
#define cf3_sdm_equations_lineuler_LinEulerUniform2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/PDE.hpp"
#include "cf3/sdm/equations/lineuler/LibLinEuler.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

/// @brief Linearized Euler 2D physics with uniform mean flow
///
/// This component assembles the terms, term-computers,
/// and configuration options
/// to solve the right-hand-side of the Euler-equations in
/// the form  dQ/dt = RHS
/// @author Willem Deconinck
class sdm_equations_lineuler_API LinEulerUniform2D : public solver::PDE {

public: // functions

  /// Contructor
  /// @param name of the component
  LinEulerUniform2D ( const std::string& name );

  /// Virtual destructor
  virtual ~LinEulerUniform2D();

  /// Get the class name
  static std::string type_name () { return "LinEulerUniform2D"; }

};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_lineuler_LinEulerUniform2D_hpp
