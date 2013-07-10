// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/euler/BCPressureOutlet.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCPressureOutlet2D,solver::BC,LibEuler> BCPressureOutlet2D_builder;

/////////////////////////////////////////////////////////////////////////////

BCPressureOutlet2D::BCPressureOutlet2D(const std::string& name) : 
  dcm::core::BC<2,4>(name),
  m_gamma(1.4)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");

  m_p = 101300.;
  options().add("p",common::to_str(m_p))
      .description("pressure as a function of x,y")
      .attach_trigger( boost::bind( &BCPressureOutlet2D::config_p, this) )
      .mark_basic();
  config_p();

}

/////////////////////////////////////////////////////////////////////////////

void BCPressureOutlet2D::config_p()
{ 
  m_function_p.parse(options().option("p").value_str(),"x,y");
}

/////////////////////////////////////////////////////////////////////////////

void BCPressureOutlet2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                                    const ColVector_NDIM& coords,
                                                    const ColVector_NDIM& face_normal,
                                                    RowVector_NEQS& boundary_solution )
{
  m_function_p.evaluate(coords,m_p);
  m_rho_inner  = inner_solution[Rho];
  m_u_inner    = inner_solution[RhoUx]/m_rho_inner;
  m_v_inner    = inner_solution[RhoUy]/m_rho_inner;
  m_uuvv_inner = m_u_inner*m_u_inner + m_v_inner*m_v_inner;

  boundary_solution[Rho  ]=inner_solution[Rho];
  boundary_solution[RhoUx]=inner_solution[RhoUx];
  boundary_solution[RhoUy]=inner_solution[RhoUy];
  boundary_solution[RhoE ]=m_p/(m_gamma-1.) + 0.5 * m_rho_inner * m_uuvv_inner;
}

/////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
