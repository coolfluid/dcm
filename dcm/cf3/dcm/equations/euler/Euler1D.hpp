// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_euler_Euler1D_hpp
#define cf3_dcm_equations_euler_Euler1D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/PDE.hpp"
#include "cf3/dcm/equations/euler/LibEuler.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

/// @brief Euler 1D physics
///
/// This component assembles the terms, term-computers,
/// and configuration options
/// to solve the right-hand-side of the Euler-equations in
/// the form  dQ/dt = RHS
/// @author Willem Deconinck
class dcm_equations_euler_API Euler1D : public solver::PDE {

public: // functions

  /// Contructor
  /// @param name of the component
  Euler1D ( const std::string& name );

  /// Virtual destructor
  virtual ~Euler1D();

  /// Get the class name
  static std::string type_name () { return "Euler1D"; }

  virtual std::string solution_variables() const;

};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_euler_Euler1D_hpp
