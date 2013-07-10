// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/dcm/equations/euler/BCPressureVelocityInlet.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

//////////////////////////////////////////////////////////////////////////////

  common::ComponentBuilder<BCPressureVelocityInlet2D,solver::BC,LibEuler> BCPressureVelocityInlet2d_builder;

/////////////////////////////////////////////////////////////////////////////

BCPressureVelocityInlet2D::BCPressureVelocityInlet2D(const std::string& name) : 
  dcm::core::BC<2,4>(name),
  m_gamma(1.4),
  m_R(287.05)
{
  m_function_p.parse("100000","x,y");  // 25 degrees Celcius
  m_function_u.parse("0","x,y");  // 1 bar


  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("R",m_R).link_to(&m_gamma)
      .mark_basic()
      .description("Gas constant");


  options().add("p",m_function_p.function()).description("Pressure")
      .attach_trigger( boost::bind( &BCPressureVelocityInlet2D::config_p, this) )
      .mark_basic();
  options().add("u",m_function_u.function()).description("Velocity magnitude")
      .attach_trigger( boost::bind( &BCPressureVelocityInlet2D::config_u, this) )
      .mark_basic();
}

/////////////////////////////////////////////////////////////////////////////

void BCPressureVelocityInlet2D::config_p()    { m_function_p.parse(options().option("p").value_str()); }
void BCPressureVelocityInlet2D::config_u()    { m_function_u.parse(options().option("u").value_str()); }

/////////////////////////////////////////////////////////////////////////////

void BCPressureVelocityInlet2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                                           const ColVector_NDIM& coords,
                                                           const ColVector_NDIM& face_normal,
                                                           RowVector_NEQS& boundary_solution )
{  
  using math::Functions::sign;

  // IMPOSE pressure and velocity:
  // -----------------------------
  // Evaluate analytical functions
  m_function_p.evaluate(coords,m_p);
  m_function_u.evaluate(coords,m_U[XX]);
  m_U[YY] = 0.;

  const Real gm1 = (m_gamma-1.);

  // USE inner density:
  // --------------------
  m_rho = inner_solution[Rho];
  
  
  m_uuvv = m_U[XX]*m_U[XX]+m_U[YY]*m_U[YY];
  m_rhoE = m_p/gm1 + 0.5*m_rho*m_uuvv;

  boundary_solution[Rho  ]=m_rho;
  boundary_solution[RhoUx]=m_rho*m_U[XX];
  boundary_solution[RhoUy]=m_rho*m_U[YY];
  boundary_solution[RhoE ]=m_p/gm1 + 0.5*m_rho*m_uuvv;;
}

/////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
