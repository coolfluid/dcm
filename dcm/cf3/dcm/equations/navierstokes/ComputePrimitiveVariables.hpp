// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_navierstokes_ComputePrimitiveVariables_hpp
#define cf3_dcm_equations_navierstokes_ComputePrimitiveVariables_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/equations/navierstokes/LibNavierStokes.hpp"
#include "cf3/common/Action.hpp"

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
namespace cf3 {
  namespace mesh {
    class Field;
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

/// @brief Compute Primitive Variables from Solution
///
/// This component computes the primitive variables from
/// the conservative variables
///
///     primitive variables:   p, U, T
/// 
///     conservative variables rho, rho*U, rho*E
///
/// @author Willem Deconinck

class dcm_equations_navierstokes_API ComputePrimitiveVariables : public common::Action {

public: // functions

  /// Contructor
  /// @param name of the component
  ComputePrimitiveVariables ( const std::string& name );

  /// Virtual destructor
  virtual ~ComputePrimitiveVariables();

  /// Get the class name
  static std::string type_name () { return "ComputePrimitiveVariables"; }

  virtual void execute();

private:
  
  Handle<mesh::Field> m_solution;
  Handle<mesh::Field> m_pressure;
  Handle<mesh::Field> m_velocity;
  Handle<mesh::Field> m_temperature;
  
};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_navierstokes_ComputePrimitiveVariables_hpp
