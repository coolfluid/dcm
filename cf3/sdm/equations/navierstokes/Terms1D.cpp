// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/bind.hpp>
#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Reconstructions.hpp"
#include "cf3/solver/TermComputer.hpp"
#include "cf3/sdm/core/CombinedTermComputer.hpp"
#include "cf3/sdm/equations/navierstokes/Terms1D.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

common::ComponentBuilder<Terms1D,solver::Term,LibNavierStokes> Terms1D_builder;
common::ComponentBuilder<core::CombinedTermComputer<Terms1D>,solver::TermComputer,LibNavierStokes> Terms1DComputer_builder;

////////////////////////////////////////////////////////////////////////////////

Terms1D::Terms1D( const std::string& name ) :
  solver::TermBase<1,3,3,2>(name),
  m_riemann_solver_type("Roe")
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("R",m_R).link_to(&m_R)
      .description("Gas constant")
      .mark_basic();
  options().add("riemann_solver",m_riemann_solver_type).link_to(&m_riemann_solver_type)
      .description("Riemann solver")
      .attach_trigger( boost::bind( &Terms1D::config_riemann_solver, this) );
  config_riemann_solver();
  options().add("k",m_k).link_to(&m_k)
      .description("Heat conduction")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Dynamic viscosity")
      .mark_basic();
}

void Terms1D::get_variables( const mesh::Space& space,
                             const Uint elem_idx,
                             const ColVector_NDIM& coords,
                             const mesh::ReconstructPoint& interpolation,
                             const std::vector<mesh::ReconstructPoint>& gradient,
                             const Matrix_NDIMxNDIM& jacobian_inverse,
                             RowVector_NVAR& vars,
                             RowVector_NGRAD& gradvars,
                             Matrix_NDIMxNGRAD& gradvars_grad )
{
  mesh::Connectivity::ConstRow nodes = space.connectivity()[elem_idx];
  vars.setZero();
  boost_foreach( Uint n, interpolation.used_points() )
  {
    const Real C = interpolation.coeff(n);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += C * solution()->array()[nodes[n]][eq];
    }
  }

  Real rho = vars[0];
  Real u = vars[1]/rho;
  Real E = vars[2]/rho;
  Real p = (m_gamma-1.)*rho*(E - 0.5*(u*u));
  Real T = p/(rho*m_R);

  gradvars[0] = u;
  gradvars[1] = T;

//  std::cout << "\nx = " << coords[XX] << std::endl;
//  std::cout << "rho  = " << rho << std::endl;
//  std::cout << "u    = " << u << std::endl;
//  std::cout << "p    = " << p << std::endl;
//  std::cout << "T    = " << T << std::endl;

  gradvars_grad.setZero();
  for (Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real C = gradient[d].coeff(n);

      Real rho = solution()->array()[nodes[n]][0];
      Real u = solution()->array()[nodes[n]][1]/rho;
      Real u2 = u*u;
      Real E = solution()->array()[nodes[n]][2]/rho;
      Real p = (m_gamma-1.)*rho*(E - 0.5*u2);
      Real T = p/(rho*m_R);

      gradvars_grad(d,0) += C * u;
      gradvars_grad(d,1) += C * T;
    }
  }
  gradvars_grad = jacobian_inverse*gradvars_grad;
}

void Terms1D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
  phys_data.mu = m_mu;
  phys_data.k = m_k;
  phys_data.Cp = m_gamma*m_R/(m_gamma-1.);
}

void Terms1D::compute_phys_data( const ColVector_NDIM& coords,
                                 const RowVector_NVAR& vars,
                                 const RowVector_NGRAD& gradvars,
                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                 DATA& phys_data )
{
  phys_data.compute_from_conservative(vars);
  phys_data.grad_u = gradvars_grad.col(0);
  phys_data.grad_T = gradvars_grad.col(1);
}

void Terms1D::compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed)
{
  physics::navierstokes::navierstokes1d::compute_diffusive_flux(p,normal,flux,wave_speed);
  std::cout << "Fd = " << flux << std::endl;

//  const Real& nx = normal[XX];
//  const Real& ny = normal[YY];

//  const Real& rho  = p.rho;
//  const Real& rhou = p.rho*p.U[XX];
//  const Real& rhov = p.rho*p.V[XX];
//  const Real& rhoE = data.solution[3];

//  const RealVectorNDIM& grad_rho  = data.solution_gradient.col(0);
//  const RealVectorNDIM& grad_rhou = data.solution_gradient.col(1);
//  const RealVectorNDIM& grad_rhov = data.solution_gradient.col(2);
//  const RealVectorNDIM& grad_rhoE = data.solution_gradient.col(3);

//  rho2 = rho*rho;
//  rho3 = rho2*rho;
//  grad_u = (rho*grad_rhou-rhou*grad_rho)/rho2;
//  grad_v = (rho*grad_rhov-rhov*grad_rho)/rho2;
//  grad_T = m_gamma_minus_1/(m_R*rho3) * (grad_rho*(-rho*rhoE+rhou*rhou+rhov*rhov)
//                                         + rho*(rho*grad_rhoE-rhou*grad_rhou-rhov*grad_rhov));
//  two_third_divergence_U = 2./3.*(grad_u[XX] + grad_v[YY]);

//  // Viscous stress tensor
//  tau_xx = m_mu*(2.*grad_u[XX] - two_third_divergence_U);
//  tau_yy = m_mu*(2.*grad_v[YY] - two_third_divergence_U);
//  tau_xy = m_mu*(grad_u[YY] + grad_v[XX]);

//  // Heat flux
//  heat_flux = -m_k*(grad_T[XX]*nx + grad_T[YY]*ny);

//  flux[0] = 0.;
//  flux[1] = tau_xx*nx + tau_xy*ny;
//  flux[2] = tau_xy*nx + tau_yy*ny;
//  flux[3] = (tau_xx*rhou + tau_xy*rhov)/rho*nx + (tau_xy*rhou + tau_yy*rhov)/rho*ny - heat_flux;

//  // maximum of kinematic viscosity nu and thermal diffusivity alpha
//  wave_speed = std::max(m_mu/rho, m_k/(rho*m_Cp));
}

void Terms1D::config_riemann_solver()
{
  std::string builder_name = m_riemann_solver_type;
  if (! boost::algorithm::starts_with(builder_name,"cf3." ) )
    builder_name = "cf3.sdm.equations.navierstokes."+builder_name+"1D";
  CFinfo << "buildername = " << builder_name << CFendl;
  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<Terms1D> >();
}
////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3
