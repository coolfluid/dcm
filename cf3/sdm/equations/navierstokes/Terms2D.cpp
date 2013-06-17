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

common::ComponentBuilder<Convection2D,solver::Term,LibNavierStokes> Convection2D_builder;
common::ComponentBuilder<core::CombinedTermComputer<Convection2D>,solver::TermComputer,LibNavierStokes> Convection2DComputer_builder;

common::ComponentBuilder<Diffusion2D,solver::Term,LibNavierStokes> Diffusion2D_builder;
common::ComponentBuilder<core::CombinedTermComputer<Diffusion2D>,solver::TermComputer,LibNavierStokes> Diffusion2DComputer_builder;

////////////////////////////////////////////////////////////////////////////////

Terms2D::Terms2D( const std::string& name ) :
  solver::TermBase<2,4,4,3>(name),
  m_riemann_solver_type("Roe"),
  m_gamma(1.4),
  m_R(287.05),
  m_kappa(2.601e-2),
  m_mu(1.806e-5)
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
  options().add("kappa",m_kappa).link_to(&m_kappa)
      .description("Thermal conduction")
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
                             const Matrix_NDIMxNDIM& jacobian,
                             const Matrix_NDIMxNDIM& jacobian_inverse,
                             const Real& jacobian_determinant,
                             RowVector_NVAR& vars,
                             RowVector_NGRAD& gradvars,
                             Matrix_NDIMxNGRAD& gradvars_grad )
{
  const mesh::Field& solution_field = *solution();
  mesh::Connectivity::ConstRow nodes = space.connectivity()[elem_idx];
  vars.setZero();

  boost_foreach( Uint n, interpolation.used_points() )
  {
    const Uint pt = nodes[n];
    const Real L = interpolation.coeff(n);

    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += L * solution_field[pt][eq];
    }
  }

  // Compute variables in point used for gradient
  _rho = vars[0];
  _u   = vars[1]/_rho;
  _v   = vars[2]/_rho;
  _E   = vars[3]/_rho;
  _U2  = _u*_u+_v*_v;
  _p   = (m_gamma-1.)*_rho*(_E - 0.5*_U2);
  _T   = _p/(_rho*m_R);

  gradvars[0] = _u;
  gradvars[1] = _v;
  gradvars[2] = _T;

  gradvars_grad.setZero();
  for (Uint d=0; d<NDIM; ++d)
  {
    const mesh::ReconstructPoint& derivative = gradient[d];
    boost_foreach( Uint n, derivative.used_points() )
    {
      const Uint pt = nodes[n];
      const Real D = derivative.coeff(n);

      _rho = solution_field[pt][0];
      _u   = solution_field[pt][1]/_rho;
      _v   = solution_field[pt][2]/_rho;
      _E   = solution_field[pt][3]/_rho;
      _U2  = _u*_u+_v*_v;
      _p   = (m_gamma-1.)*_rho*(_E - 0.5*_U2);
      _T   = _p/(_rho*m_R);

      gradvars_grad(d,0) += D * _u;
      gradvars_grad(d,1) += D * _v;
      gradvars_grad(d,2) += D * _T;
    }
  }
  gradvars_grad = jacobian_inverse * gradvars_grad;
}

void Terms2D::get_bdry_variables( const mesh::Space& space,
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
    const Uint pt = nodes[n];
    const Real L = interpolation.coeff(n);

    for (Uint eq=0; eq<NEQS; ++eq)
      vars[eq] += L * bdry_solution_field[pt][eq];
  }

  // Compute variables in point used for gradient
  _rho = vars[0];
  _u   = vars[1]/_rho;
  _v   = vars[2]/_rho;
  _E   = vars[3]/_rho;
  _U2  = _u*_u+_v*_v;
  _p   = (m_gamma-1.)*_rho*(_E - 0.5*_U2);
  _T   = _p/(_rho*m_R);

  for (Uint d=0; d<NDIM; ++d)
  {
    const mesh::ReconstructPoint& interpolate_derivative = gradient[d];
    boost_foreach( Uint n, interpolate_derivative.used_points() )
    {
      const Uint pt = nodes[n];
      const Real L = interpolate_derivative.coeff(n);
      _grad_rho[d]  = L * bdry_solution_gradient_field[pt][0+d*NEQS];
      _grad_rhou[d] = L * bdry_solution_gradient_field[pt][1+d*NEQS];
      _grad_rhov[d] = L * bdry_solution_gradient_field[pt][2+d*NEQS];
      _grad_rhoE[d] = L * bdry_solution_gradient_field[pt][3+d*NEQS];
    }
  }

  _grad_u = 1./_rho * ( _grad_rhou - _u*_grad_rho );
  _grad_v = 1./_rho * ( _grad_rhov - _v*_grad_rho );
  _grad_T = ( (_grad_rhoE - _grad_u*_u - _grad_v*_v)*(m_gamma-1)*_rho -_p*_grad_rho )/(m_R*_rho*_rho);

  gradvars[0] = _u;
  gradvars[1] = _v;
  gradvars[2] = _T;
  gradvars_grad.col(0) = _grad_u;
  gradvars_grad.col(1) = _grad_v;
  gradvars_grad.col(2) = _grad_T;
}


void Terms2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
  phys_data.mu = m_mu;
  phys_data.kappa = m_kappa;
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
  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<Terms2D> >();
}
////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3
