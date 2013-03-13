// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_navierstokes_InitBlasius_hpp
#define cf3_sdm_equations_navierstokes_InitBlasius_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Action.hpp"
#include "cf3/sdm/equations/navierstokes/LibNavierStokes.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace mesh { class Field; }
namespace sdm { 
namespace equations {
namespace navierstokes {

//////////////////////////////////////////////////////////////////////////////

/// This configurable function creates an initial condition for a
/// Navier-Stokes Blasius profile over flat plate
/// @author Willem Deconinck
class  InitBlasius : public common::Action
{
public: // functions
  
  /// constructor
  InitBlasius( const std::string& name );
  
  virtual ~InitBlasius() {}
  /// Gets the Class name
  static std::string type_name() { return "InitBlasius"; }

  virtual void execute();

  void compute_blasius( const std::vector<Real> coords, std::vector<Real>& cons );
  
private: // data

  Handle<mesh::Field> m_field;
  Real m_M_inf;
  Real m_T_inf;
  Real m_p_inf;
  Real m_L;
  Real m_Re;
  Real m_Pr;
  Real m_R;
  Real m_gamma;
  
  
  Real m_c_inf;
  Real m_u_inf;
  Real m_rho_inf;
  Real m_k;
  Real m_mu;
  Real m_nu;
  Real m_Cp;
  Real m_rhoE_inf;

  
}; // end InitBlasius


////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_navierstokes_InitBlasius_hpp
