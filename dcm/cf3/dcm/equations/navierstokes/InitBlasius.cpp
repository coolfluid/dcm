// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

//#include <boost/math/special_functions/bessel.hpp>

#include "cf3/common/Log.hpp"
#include "cf3/common/Builder.hpp"

#include "cf3/common/FindComponents.hpp"
#include "cf3/common/Foreach.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionT.hpp"
#include "cf3/common/OptionComponent.hpp"

#include "cf3/mesh/Elements.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Space.hpp"

#include "cf3/dcm/equations/navierstokes/InitBlasius.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

  using namespace common;

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < InitBlasius, common::Action, LibNavierStokes> InitBlasius_Builder;

//////////////////////////////////////////////////////////////////////////////

InitBlasius::InitBlasius( const std::string& name )
  : common::Action(name)
{

  properties()["brief"] = std::string("Initialize a field with a constant value");
  std::string desc;
  desc = "  Usage: InitBlasius constant \n";
  properties()["description"] = desc;

  options().add("field", m_field)
      .description("Field to initialize")
      .pretty_name("Field")
      .link_to(&m_field)
      .mark_basic();

  // AIR
  m_R = 287.05;
  m_gamma = 1.4;
  m_Pr = 0.72;
  
  m_T_inf = 0.;
  m_M_inf = 0.;
  m_p_inf = 0.;
  
  m_Re = 0.;
  m_L = 0.;

  options().add("R",m_R).link_to(&m_R).mark_basic();
  options().add("gamma",m_gamma).link_to(&m_gamma).mark_basic();
  options().add("Re",m_Re).link_to(&m_Re).mark_basic();
  options().add("Pr",m_Pr).link_to(&m_Pr).mark_basic();
  options().add("L",m_L).link_to(&m_L).mark_basic();
  options().add("M",m_M_inf).link_to(&m_M_inf).mark_basic();
  options().add("T",m_T_inf).link_to(&m_T_inf).mark_basic();
  options().add("p",m_p_inf).link_to(&m_p_inf).mark_basic();
}

////////////////////////////////////////////////////////////////////////////////

void InitBlasius::execute()
{
  if (is_null(m_field)) throw SetupError( FromHere(), "field not configured" );
  if (m_L == 0.) throw SetupError( FromHere(), "L (plate length) not configured" );
  if (m_Re == 0.) throw SetupError( FromHere(), "Re (Reynolds) not configured" );
  if (m_T_inf == 0.) throw SetupError( FromHere(), "T (temperature) not configured" );
  if (m_p_inf == 0.) throw SetupError( FromHere(), "p (pressure) not configured" );
  if (m_M_inf == 0.) throw SetupError( FromHere(), "M (Mach) not configured" );

  m_rho_inf = m_p_inf/(m_R*m_T_inf);
  m_c_inf = std::sqrt(m_gamma*m_R*m_T_inf);
  m_u_inf = m_M_inf * m_c_inf;
  m_nu = m_u_inf * m_L / m_Re;
  m_mu = m_rho_inf * m_nu;
  m_Cp = m_gamma*m_R/(m_gamma-1.);
  m_k = m_Cp * m_mu / m_Pr;
  m_rhoE_inf = m_p_inf/(m_gamma-1) + 0.5*m_rho_inf*(m_u_inf*m_u_inf);
      
  std::vector<Real> coord(DIM_2D);
  std::vector<Real> cons(4);

  cf3_assert(m_field->coordinates().row_size()>=DIM_2D);
  cf3_assert(m_field->row_size()==4);
  for (Uint i=0; i<m_field->size(); ++i)
  {
    cf3_assert(i<m_field->coordinates().size());
    coord[XX] = m_field->coordinates()[i][XX];
    coord[YY] = m_field->coordinates()[i][YY];

    compute_blasius(coord, cons);

    m_field->array()[i][0] = cons[0];
    m_field->array()[i][1] = cons[1];
    m_field->array()[i][2] = cons[2];
    m_field->array()[i][3] = cons[3];
  }

}

//////////////////////////////////////////////////////////////////////////////


void InitBlasius::compute_blasius( const std::vector<Real> coords, std::vector<Real>& cons )
{
  const Real x = coords[XX];
  const Real y = coords[YY];

  // upstream conditions before flat plate, including the leading edge
  if ( x <= 0. || y > 0.5*m_L)
  {
    cons[0] = m_rho_inf;
    cons[1] = m_rho_inf * m_u_inf;
    cons[2] = 0.;
    cons[3] = m_rhoE_inf; 
    return;
  }

  Real f, fo, fp, fpp;
  Real eta;
  Real n, dn, k1, k2, k3, k4;

  eta=0.;  f=0.; fo=0.; fp=0.; fpp=0.33206;
  n=0.; dn=0.00005;

  // Determine the dimensionless similarity coordinate, eta:
  eta = y*std::sqrt(m_u_inf/(x*m_mu/m_rho_inf));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is
  // linear in this region)

  if (eta > 100.0) {
    cons[0] = m_rho_inf;
    cons[1] = m_rho_inf * m_u_inf;
    cons[2] = 0.;
    cons[3] = m_rhoE_inf; 
    return;
  } else if (eta > 8.4) {
    fp = 1.; fpp = 0.; f = 6.67923 + fp*(eta - 8.4);
    const Real u = fp*m_u_inf;
    const Real v = std::sqrt(0.5*m_nu*m_u_inf/x)*(eta*fp-f);
    cons[0] = m_rho_inf;
    cons[1] = m_rho_inf * u;
    cons[2] = m_rho_inf * v;
    cons[3] = m_p_inf/(m_gamma-1) + 0.5*m_rho_inf*(u*u+v*v);
    return;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method
  while (n < eta) {

    // Store the solution at the start of the iteration
    fo = f;

    // Increment f:
    k1 = dn*fp;
    k2 = dn*(fp + k1/2.0);
    k3 = dn*(fp + k2/2.0);
    k4 = dn*(fp + k3);
    f += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fp:
    k1 = dn*fpp;
    k2 = dn*(fpp + k1/2.0);
    k3 = dn*(fpp + k2/2.0);
    k4 = dn*(fpp + k3);
    fp += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fpp:
    k1 = -dn*fo*fpp/2.0;
    k2 = -dn*(fo + dn/2.0)*(fpp + k1/2.0)/2.0;
    k3 = -dn*(fo + dn/2.0)*(fpp + k2/2.0)/2.0;
    k4 = -dn*(fo + dn)*(fpp + k3)/2.0;
    fpp += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Determine the increment dn:
    if (n + dn > eta) dn = eta - n;

    // Increment n:
    n += dn;
  }

  // Compute the velocity vector at point X
  const Real u = fp*m_u_inf;
  const Real v = std::sqrt(0.5*m_nu*m_u_inf/x)*(eta*fp-f);
  cons[0] = m_rho_inf;
  cons[1] = m_rho_inf * u;
  cons[2] = m_rho_inf * v;
  cons[3] = m_p_inf/(m_gamma-1) + 0.5*m_rho_inf*(u*u+v*v);
}


} // navierstokes
} // equations
} // dcm
} // cf3
