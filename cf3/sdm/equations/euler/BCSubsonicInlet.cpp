// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/sdm/equations/euler/BCSubsonicInlet.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

//////////////////////////////////////////////////////////////////////////////

  common::ComponentBuilder<BCSubsonicInlet2D,solver::BC,LibEuler> bcsubsonicinletTtPtAlpha2d_builder;

/////////////////////////////////////////////////////////////////////////////

BCSubsonicInlet2D::BCSubsonicInlet2D(const std::string& name) : 
  sdm::core::BC<2,4>(name),
  m_gamma(1.4),
  m_R(287.05)
{
  m_function_Tt.parse("298.15","x,y");  // 25 degrees Celcius
  m_function_Pt.parse("100000","x,y");  // 1 bar
  m_function_alpha.parse("0","x,y");    // 0 radians


  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("R",m_R).link_to(&m_gamma)
      .mark_basic()
      .description("Gas constant");


  options().add("Tt",m_function_Tt.function()).description("Total Temperature")
      .attach_trigger( boost::bind( &BCSubsonicInlet2D::config_Tt, this) )
      .mark_basic();
  options().add("Pt",m_function_Pt.function()).description("Total Pressure")
      .attach_trigger( boost::bind( &BCSubsonicInlet2D::config_Pt, this) )
      .mark_basic();
  options().add("alpha",m_function_alpha.function()).description("flow angle in rad")
      .attach_trigger( boost::bind( &BCSubsonicInlet2D::config_alpha, this) )
      .mark_basic();
}

/////////////////////////////////////////////////////////////////////////////

void BCSubsonicInlet2D::config_Tt()    { m_function_Tt   .parse(options().option("Tt").value_str()); }
void BCSubsonicInlet2D::config_Pt()    { m_function_Pt   .parse(options().option("Pt").value_str()); }
void BCSubsonicInlet2D::config_alpha() { m_function_alpha.parse(options().option("alpha").value_str()); }

/////////////////////////////////////////////////////////////////////////////

void BCSubsonicInlet2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                                            const ColVector_NDIM& coords,
                                                            const ColVector_NDIM& face_normal,
                                                            RowVector_NEQS& boundary_solution )
{  
  using math::Functions::sign;

  // Evaluate analytical functions
  m_function_Tt.evaluate(coords,m_Tt);
  m_function_Pt.evaluate(coords,m_Pt);
  m_function_alpha.evaluate(coords,m_alpha);

  const Real gm1 = (m_gamma-1.);
  // Compute inner cell data
  m_rho_inner       = inner_solution[Rho];
  m_uuvv_inner      = (inner_solution[RhoUx]*inner_solution[RhoUx] + inner_solution[RhoUy]*inner_solution[RhoUy])/(m_rho_inner*m_rho_inner);
  m_rhoE_inner      = inner_solution[RhoE];
  m_p_inner         = gm1*(m_rhoE_inner - 0.5 * m_rho_inner * m_uuvv_inner);
  m_T_inner         = m_p_inner / (m_R*m_rho_inner);
  m_M2_inner        = m_uuvv_inner/(m_gamma*m_R*m_T_inner);
  m_coeff_inner     = 1. + 0.5*gm1*m_M2_inner;
  m_pow_coeff_inner = std::pow(m_coeff_inner,m_gamma/gm1);
  //m_Tt_inner    = m_T_inner*m_coeff_inner;
  //m_Pt_inner    = m_p_inner*m_pow_coeff_inner;

  // Compute values to impose on boundary
  m_M = sqrt(m_M2_inner);
  m_tan_alpha=std::tan(m_alpha);
  m_T = m_Tt/m_coeff_inner;
  m_p = m_Pt/m_pow_coeff_inner;
  m_rho = m_p/(m_R*m_T);
  m_U[XX] = sign(std::cos(m_alpha)) * m_M*std::sqrt(m_gamma*m_R*m_T/(1.+m_tan_alpha*m_tan_alpha));
  m_U[YY] = m_tan_alpha*m_U[XX];
  m_uuvv = m_U[XX]*m_U[XX]+m_U[YY]*m_U[YY];
  m_rhoE = m_p/gm1 + 0.5*m_rho*m_uuvv;

  boundary_solution[Rho  ]=m_rho;
  boundary_solution[RhoUx]=m_rho*m_U[XX];
  boundary_solution[RhoUy]=m_rho*m_U[YY];
  boundary_solution[RhoE ]=m_rhoE;
}

/////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3
