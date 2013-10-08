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
#include "cf3/dcm/equations/advectiondiffusion/RightHandSide1D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

common::ComponentBuilder<RightHandSide1D,solver::Term,LibAdvectionDiffusion> RightHandSide1D_builder;

////////////////////////////////////////////////////////////////////////////////

RightHandSide1D::RightHandSide1D( const std::string& name ) :
  solver::TermBase<1,1,1,1>(name)
{
  options().add("a",m_a).link_to(&m_a)
      .description("Advection speed")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Diffusion coefficient")
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
    const Real L = interpolation.coeff(n);
    vars[0] += L * solution()->array()[nodes[n]][0];
  }
  gradvars[0] = vars[0];
  gradvars_grad.setZero();
  boost_foreach( Uint n, gradient[0].used_points() )
  {
    const Real D = gradient[0].coeff(n);
    gradvars_grad(0,0) += D * solution()->array()[nodes[n]][0];
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
    const Uint pt = nodes[n];
    const Real L = interpolation.coeff(n);
    vars[0] += L * bdry_solution_field[pt][0];
  }
  gradvars = vars;

  const mesh::ReconstructPoint& interpolate_derivative = gradient[0];
  boost_foreach( Uint n, interpolate_derivative.used_points() )
  {
    const Uint pt = nodes[n];
    const Real L = interpolate_derivative.coeff(n);
    gradvars_grad[0] = L * bdry_solution_gradient_field[pt][0];
  }
}


void RightHandSide1D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.a = m_a;
  phys_data.mu = m_mu;
}

void RightHandSide1D::compute_phys_data( const ColVector_NDIM& coords,
                                 const RowVector_NVAR& vars,
                                 const RowVector_NGRAD& gradvars,
                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                 DATA& phys_data )
{
  phys_data.q = vars[0];
  phys_data.grad_q = gradvars_grad[0];
}

void RightHandSide1D::compute_convective_flux( const DATA& p, const ColVector_NDIM& normal,
                                       RowVector_NEQS& flux, Real& wave_speed )
{
  flux[0] = p.a * normal[0] * p.q;
  wave_speed = std::abs(p.a);
}

void RightHandSide1D::compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                                    RowVector_NEQS& flux, Real& wave_speed )
{
  const Real an = left.a * normal[0];
  if (an>0) 
    flux[0] = an * left.q;
  else
    flux[0] = an * right.q;  
  wave_speed = std::abs(an);
}

void RightHandSide1D::compute_diffusive_flux( const DATA& p, const ColVector_NDIM& normal,
                                      RowVector_NEQS& flux, Real& wave_speed )
{
  flux[0] = p.mu * p.grad_q * normal[0];
  wave_speed = std::abs(p.mu);
}

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3
