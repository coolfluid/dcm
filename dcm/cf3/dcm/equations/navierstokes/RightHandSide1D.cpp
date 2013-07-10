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
#include "cf3/dcm/equations/navierstokes/RightHandSide1D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

common::ComponentBuilder<RightHandSide1D,solver::Term,LibNavierStokes> RightHandSide1D_builder;

////////////////////////////////////////////////////////////////////////////////

RightHandSide1D::RightHandSide1D( const std::string& name ) :
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
      .attach_trigger( boost::bind( &RightHandSide1D::config_riemann_solver, this) );
  config_riemann_solver();
  options().add("kappa",m_kappa).link_to(&m_kappa)
      .description("Thermal conduction")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Dynamic viscosity")
      .mark_basic();
}

void RightHandSide1D::get_variables( const mesh::Space& space,
                             const Uint elem_idx,
                             const ColVector_NDIM& coords,
                             const mesh::ReconstructPoint& interpolation,
                             const std::vector<mesh::ReconstructPoint>& gradient,
                             const Matrix_NDIMxNDIM& jacobian,
                             const Matrix_NDIMxNDIM& jacobian_inverse,
                             const Real& jacobian_determinant,
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
  gradvars_grad = jacobian_inverse * gradvars_grad;
}


void RightHandSide1D::get_bdry_variables( const mesh::Space& space,
                                  const Uint elem_idx,
                                  const ColVector_NDIM& coords,
                                  const mesh::ReconstructPoint& interpolation,
                                  const std::vector<mesh::ReconstructPoint>& gradient,
                                  const Matrix_NDIMxNDIM& jacobian,
                                  const Matrix_NDIMxNDIM& jacobian_inverse,
                                  const Real& jacobian_determinant,
                                  RowVector_NVAR& vars,
                                  RowVector_NGRAD& gradvars,
                                  Matrix_NDIMxNGRAD& gradvars_grad )
{
  const mesh::Field& bdry_solution_field = *bdry_solution();
  const mesh::Field& bdry_solution_gradient_field = *bdry_solution_gradient();
  mesh::Connectivity::ConstRow nodes = space.connectivity()[elem_idx];
  vars.setZero();
  boost_foreach( Uint n, interpolation.used_points() )
  {
    const Real C = interpolation.coeff(n);

    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += C * bdry_solution_field[nodes[n]][eq];
    }
  }
  // Compute variables in point used for gradient
  const Real rho = vars[0];
  const Real u   = vars[1]/rho;
  const Real E   = vars[2]/rho;
  const Real U2  = u*u;
  const Real p   = (m_gamma-1.)*rho*(E - 0.5*U2);
  const Real T   = p/(rho*m_R);

  gradvars[0] = u;
  gradvars[1] = T;

  ColVector_NDIM grad_rho;
  ColVector_NDIM grad_rhou;
  ColVector_NDIM grad_rhoE;
  ColVector_NDIM grad_u;
  ColVector_NDIM grad_v;
  ColVector_NDIM grad_T;

  for (Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real C = gradient[d].coeff(n);

      grad_rho[d]  = C * bdry_solution_gradient_field[nodes[n]][0+d*NEQS];
      grad_rhou[d] = C * bdry_solution_gradient_field[nodes[n]][1+d*NEQS];
      grad_rhoE[d] = C * bdry_solution_gradient_field[nodes[n]][2+d*NEQS];
    }
  }

  grad_u = 1./rho * ( grad_rhou - u*grad_rho );
  grad_T = ( (grad_rhoE - grad_u*u)*(m_gamma-1)*rho -p*grad_rho )/(m_R*rho*rho);

  gradvars_grad.col(0) = grad_u;
  gradvars_grad.col(1) = grad_T;
}



void RightHandSide1D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
  phys_data.mu = m_mu;
  phys_data.kappa = m_kappa;
  phys_data.Cp = m_gamma*m_R/(m_gamma-1.);
}

void RightHandSide1D::compute_phys_data( const ColVector_NDIM& coords,
                                 const RowVector_NVAR& vars,
                                 const RowVector_NGRAD& gradvars,
                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                 DATA& phys_data )
{
  phys_data.compute_from_conservative(vars);
  phys_data.grad_u = gradvars_grad.col(0);
  phys_data.grad_T = gradvars_grad.col(1);
}

void RightHandSide1D::compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed)
{
  physics::navierstokes::navierstokes1d::compute_diffusive_flux(p,normal,flux,wave_speed);
}

void RightHandSide1D::config_riemann_solver()
{
  std::string builder_name = m_riemann_solver_type;
  if (! boost::algorithm::starts_with(builder_name,"cf3." ) )
    builder_name = "cf3.dcm.equations.euler."+builder_name+"1D";
  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<cf3::physics::euler::euler1d::Data,NDIM,NEQS> >();
}
////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3
