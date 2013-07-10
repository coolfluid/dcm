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
#include "cf3/dcm/equations/advectiondiffusion/RightHandSide2D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

common::ComponentBuilder<RightHandSide2D,solver::Term,LibAdvectionDiffusion> RightHandSide2D_builder;

////////////////////////////////////////////////////////////////////////////////

RightHandSide2D::RightHandSide2D( const std::string& name ) :
  solver::TermBase<2,1,1,1>(name)
{
  m_a.resize(NDIM);
  options().add("a",m_a).link_to(&m_a)
      .description("Advection speed")
      .mark_basic();
  options().add("mu",m_mu).link_to(&m_mu)
      .description("Diffusion coefficient")
      .mark_basic();
}

void RightHandSide2D::get_variables( const mesh::Space& space,
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
  for (Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real C = gradient[d].coeff(n);
      gradvars_grad(d,0) += C * solution()->array()[nodes[n]][0];
    }
  }
  gradvars_grad = jacobian_inverse * gradvars_grad;
}

void RightHandSide2D::get_bdry_variables( const mesh::Space& space,
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
    vars[0] += L * bdry_solution()->array()[nodes[n]][0];
  }

  gradvars = vars;
  gradvars_grad.setZero();
  for( Uint d=0; d<NDIM; ++d)
  {
    boost_foreach( Uint n, gradient[d].used_points() )
    {
      const Real D = gradient[d].coeff(n);
      gradvars_grad[0] += D * bdry_solution()->array()[nodes[n]][0];
    }
  }
}

void RightHandSide2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.a[XX] = m_a[XX];
  phys_data.a[YY] = m_a[YY];
  phys_data.mu = m_mu;
}

void RightHandSide2D::compute_phys_data( const ColVector_NDIM& coords,
                                 const RowVector_NVAR& vars,
                                 const RowVector_NGRAD& gradvars,
                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                 DATA& phys_data )
{
  phys_data.q = vars[0];
  phys_data.grad_q = gradvars_grad;
}

void RightHandSide2D::compute_convective_flux( const DATA& p, const ColVector_NDIM& normal,
                                       RowVector_NEQS& flux, Real& wave_speed )
{
  Real an = p.a.dot(normal);
  flux[0] = an * p.q;
  wave_speed = std::abs(an);
}

void RightHandSide2D::compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                                    RowVector_NEQS& flux, Real& wave_speed )
{
  const Real an = left.a.dot(normal);
  if (an>0) 
    flux[0] = an * left.q;
  else
    flux[0] = an * right.q;  
  wave_speed = std::abs(an);
}

void RightHandSide2D::compute_diffusive_flux( const DATA& p, const ColVector_NDIM& normal,
                                      RowVector_NEQS& flux, Real& wave_speed )
{
  flux[0] = p.mu * p.grad_q.dot(normal);
  wave_speed = std::abs(p.mu);
}

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3
