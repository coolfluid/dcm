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
#include "cf3/sdm/equations/navierstokes/Terms2D.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

common::ComponentBuilder<Terms2D,solver::Term,LibNavierStokes> Terms2D_builder;
common::ComponentBuilder<core::CombinedTermComputer<Terms2D>,solver::TermComputer,LibNavierStokes> Terms2DComputer_builder;

////////////////////////////////////////////////////////////////////////////////

Terms2D::Terms2D( const std::string& name ) :
  solver::TermBase<2,4,4,3>(name),
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
      .attach_trigger( boost::bind( &Terms2D::config_riemann_solver, this) );
  config_riemann_solver();
  options().add("k",m_k).link_to(&m_k)
      .description("Heat conduction")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Dynamic viscosity")
      .mark_basic();
}

void Terms2D::get_variables( const mesh::Space& space,
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
  Real v = vars[2]/rho;
  Real E = vars[3]/rho;
  Real p = (m_gamma-1.)*rho*(E - 0.5*(u*u+v*v));
  Real T = p/(rho*m_R);

  gradvars[0] = u;
  gradvars[1] = v;
  gradvars[2] = T;

  gradvars_grad.setZero();
  for (Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real C = gradient[d].coeff(n);

      Real rho = solution()->array()[nodes[n]][0];
      Real u = solution()->array()[nodes[n]][1]/rho;
      Real v = solution()->array()[nodes[n]][2]/rho;
      Real E = solution()->array()[nodes[n]][3]/rho;
      Real U2 = u*u+v*v;
      Real p = (m_gamma-1.)*rho*(E - 0.5*U2);
      Real T = p/(rho*m_R);

      gradvars_grad(d,0) += C * u;
      gradvars_grad(d,1) += C * v;
      gradvars_grad(d,2) += C * T;
    }
  }
  gradvars_grad = jacobian_inverse*gradvars_grad;
}

void Terms2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
  phys_data.mu = m_mu;
  phys_data.k = m_k;
  phys_data.Cp = m_gamma*m_R/(m_gamma-1.);
}

void Terms2D::compute_phys_data( const ColVector_NDIM& coords,
                                 const RowVector_NVAR& vars,
                                 const RowVector_NGRAD& gradvars,
                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                 DATA& phys_data )
{
  phys_data.compute_from_conservative(vars);
  phys_data.grad_u = gradvars_grad.col(0);
  phys_data.grad_v = gradvars_grad.col(1);
  phys_data.grad_T = gradvars_grad.col(2);
}

void Terms2D::compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed)
{
  physics::navierstokes::navierstokes2d::compute_diffusive_flux(p,normal,flux,wave_speed);
}

void Terms2D::config_riemann_solver()
{
  std::string builder_name = m_riemann_solver_type;
  if (! boost::algorithm::starts_with(builder_name,"cf3." ) )
    builder_name = "cf3.sdm.equations.navierstokes."+builder_name+"2D";
  CFinfo << "buildername = " << builder_name << CFendl;
  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<Terms2D> >();
}
////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3
