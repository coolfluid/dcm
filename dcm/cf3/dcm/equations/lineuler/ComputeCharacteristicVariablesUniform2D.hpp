// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_ComputeCharacteristicVariablesUniform2D_hpp
#define cf3_dcm_equations_lineuler_ComputeCharacteristicVariablesUniform2D_hpp

#include "cf3/common/Action.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"

/////////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace mesh   { class Field; }
namespace dcm {
namespace equations {
namespace lineuler {

class dcm_equations_lineuler_API ComputeCharacteristicVariablesUniform2D : public common::Action
{
public: // typedefs

  /// pointers
  typedef boost::shared_ptr<ComputeCharacteristicVariablesUniform2D> Ptr;
  typedef boost::shared_ptr<ComputeCharacteristicVariablesUniform2D const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  ComputeCharacteristicVariablesUniform2D ( const std::string& name );

  /// Virtual destructor
  virtual ~ComputeCharacteristicVariablesUniform2D() {};

  /// Get the class name
  static std::string type_name () { return "ComputeCharacteristicVariablesUniform2D"; }

  /// execute the action
  virtual void execute ();

private: // data

  Handle<mesh::Field> m_cons;
  Handle<mesh::Field> m_char;
  Handle<mesh::Field> m_gradn_char;
  
  Real m_c0;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

/////////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_ComputeCharacteristicVariables_hpp
