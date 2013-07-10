// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_les_LES2D_hpp
#define cf3_dcm_equations_les_LES2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/PDE.hpp"
#include "cf3/dcm/equations/les/LibLES.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
  namespace mesh {
    class Field;
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

/// @brief LES 2D physics
///
/// This component assembles the terms, term-computers,
/// and configuration options
/// to solve the right-hand-side of the Large Eddy Simulation-equations in
/// the form  dQ/dt = RHS
/// @author Willem Deconinck
class dcm_equations_les_API LES2D : public solver::PDE {

public: // functions

  /// Contructor
  /// @param name of the component
  LES2D ( const std::string& name );

  /// Virtual destructor
  virtual ~LES2D();

  /// Get the class name
  static std::string type_name () { return "LES2D"; }

  virtual std::string solution_variables() const;

private: // data
  
  Handle<mesh::Field> m_sfs_length;
    
};

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_les_LES2D_hpp
