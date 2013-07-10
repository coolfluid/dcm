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
#include "cf3/dcm/equations/les/Smagorinsky2D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

common::ComponentBuilder<Smagorinsky2D,solver::Term,LibLES> Smagorinsky2D_builder;

////////////////////////////////////////////////////////////////////////////////

Smagorinsky2D::Smagorinsky2D( const std::string& name ) :
  solver::TermBase<2,4,5,3>(name),
  m_riemann_solver_type("Roe"),
  m_PrT(0.9),
  m_Cs(0.18),
  m_Cv(0.094),
  m_gamma(1.4),
  m_kappa(2.601e-2),
  m_mu(1.806e-5),
  m_R(287.05)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("R",m_R).link_to(&m_R)
      .description("Gas constant")
      .mark_basic();
  options().add("riemann_solver",m_riemann_solver_type).link_to(&m_riemann_solver_type)
      .description("Riemann solver")
      .attach_trigger( boost::bind( &Smagorinsky2D::config_riemann_solver, this) );
  config_riemann_solver();
  options().add("kappa",m_kappa).link_to(&m_kappa)
      .description("Heat conduction")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Dynamic viscosity")
      .mark_basic();

  // Subfilter scale parameters to tune
  options().add("PrT",m_PrT).link_to(&m_PrT)
      .description("Turbulent Prandtl number")
      .mark_basic();
  options().add("Cs",m_Cs).link_to(&m_Cs)
      .description("Smagorinsky constant")
      .mark_basic();
  options().add("Cv",m_Cv).link_to(&m_Cv)
      .description("Deardorff/Yoshizawa constant")
      .mark_basic();
}

void Smagorinsky2D::get_variables( const mesh::Space& space,
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
    const Real L = interpolation.coeff(n);

    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += L * solution()->array()[nodes[n]][eq];
    }
  }
  Real dx = std::pow(jacobian_determinant/space.shape_function().nb_nodes(), 1./static_cast<Real>(NDIM));
  vars[4] = dx;

  _rho = vars[0];
  _u = vars[1]/_rho;
  _v = vars[2]/_rho;
  _E = vars[3]/_rho;
  _p = (m_gamma-1.)*_rho*(_E - 0.5*(_u*_u+_v*_v));
  _T = _p/(_rho*m_R);

  gradvars[0] = _u;
  gradvars[1] = _v;
  gradvars[2] = _T;

  gradvars_grad.setZero();
  for (Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real D = gradient[d].coeff(n);

      _rho = solution()->array()[nodes[n]][0];
      _u = solution()->array()[nodes[n]][1]/_rho;
      _v = solution()->array()[nodes[n]][2]/_rho;
      _E = solution()->array()[nodes[n]][3]/_rho;
      _U2 = _u*_u+_v*_v;
      _p = (m_gamma-1.)*_rho*(_E - 0.5*_U2);
      _T = _p/(_rho*m_R);

      gradvars_grad(d,0) += D * _u;
      gradvars_grad(d,1) += D * _v;
      gradvars_grad(d,2) += D * _T;
    }
  }
  gradvars_grad = jacobian_inverse * gradvars_grad;
}

void Smagorinsky2D::get_bdry_variables( const mesh::Space& space,
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
  Real dx = std::pow(jacobian_determinant/space.shape_function().nb_nodes(), 1./static_cast<Real>(NDIM));
  vars[4] = dx;

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

void Smagorinsky2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
  phys_data.mu = m_mu;
  phys_data.kappa = m_kappa;
  phys_data.Cp = m_gamma*m_R/(m_gamma-1.);
  phys_data.Cs = m_Cs;
  phys_data.Cv = m_Cv;
  phys_data.PrT = m_PrT;
}

void Smagorinsky2D::compute_phys_data( const ColVector_NDIM& coords,
                                       const RowVector_NVAR& vars,
                                       const RowVector_NGRAD& gradvars,
                                       const Matrix_NDIMxNGRAD& gradvars_grad,
                                       DATA& p )
{
  RowVector_NEQS sol; sol << vars[0], vars[1], vars[2], vars[3];
  p.compute_from_conservative(sol);
  p.grad_u = gradvars_grad.col(0);
  p.grad_v = gradvars_grad.col(1);
  p.grad_T = gradvars_grad.col(2);

  return;
  p.Delta  = vars[4];

  // Strain rate tensor
  Real Sxx = p.grad_u[XX];
  Real Sxy = 0.5*(p.grad_u[YY]+p.grad_v[XX]);
  Real Syy = p.grad_v[YY];

  // Absolute strain rate tensor
  Real SijSij = Sxx*Sxx + Syy*Syy + 2*Sxy*Sxy;
  Real absS = std::sqrt( 2.* SijSij);

  // Eddy viscosity computed by Smagorinsky model
  p.nuT = std::pow(p.Delta*p.Cs, 2)*absS;
  p.muT = p.rho*p.nuT;

  // SFS thermal conductivity
  p.kappaT = p.muT/p.PrT;

  // Subfilterscale kinetic energy model by Yoshizawa, Deardorff
  p.k_sfs = std::pow( p.nuT/(p.Cv*p.Delta) , 2);
}

void Smagorinsky2D::compute_diffusive_flux( const DATA &p, const ColVector_NDIM &normal,
                                            RowVector_NEQS &flux, Real &wave_speed )
{
  physics::navierstokes::navierstokes2d::compute_diffusive_flux(p,normal,flux,wave_speed);
  return;
  const Real& nx = normal[XX];
  const Real& ny = normal[YY];

  Real two_third_divergence_U = 2./3.*(p.grad_u[XX] + p.grad_v[YY]);

  Real two_third_rho_ksfs = 2./3.*p.rho*p.k_sfs;

  // Viscous stress tensor + SFS stress tensor
  // tau_ij = (mu + muT) ( du_i/dx_j + du_j/dx_i - delta_ij 2/3 div(u) ) - delta_ij 2/3 rho k_sfs
  Real mu = p.mu + p.muT;
  Real tau_xx = mu*(2.*p.grad_u[XX] - two_third_divergence_U) - two_third_rho_ksfs;
  Real tau_yy = mu*(2.*p.grad_v[YY] - two_third_divergence_U) - two_third_rho_ksfs;
  Real tau_xy = mu*(p.grad_u[YY] + p.grad_v[XX]);

  // Heat flux + SFS heat flus
  Real kappa = p.kappa + p.kappaT;
  Real heat_flux = -kappa*(p.grad_T[XX]*nx + p.grad_T[YY]*ny);

  flux[0] = 0.;
  flux[1] = tau_xx*nx + tau_xy*ny;
  flux[2] = tau_xy*nx + tau_yy*ny;
  flux[3] = (tau_xx*p.U[XX] + tau_xy*p.U[YY])*nx + (tau_xy*p.U[XX] + tau_yy*p.U[YY])*ny - heat_flux;

  wave_speed = std::max(mu/p.rho, kappa/(p.rho*p.Cp));

}


void Smagorinsky2D::config_riemann_solver()
{
  std::string builder_name = m_riemann_solver_type;
  if (! boost::algorithm::starts_with(builder_name,"cf3." ) )
    builder_name = "cf3.dcm.equations.euler."+builder_name+"2D";
  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<cf3::physics::euler::euler2d::Data,NDIM,NEQS> >();
}
////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3
