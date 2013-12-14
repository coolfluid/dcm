// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Reconstructions.hpp"
#include "cf3/dcm/equations/lineuler/RightHandSide2D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

common::ComponentBuilder<RightHandSide2D,solver::Term,LibLinEuler> RightHandSide2D_builder;

////////////////////////////////////////////////////////////////////////////////

RightHandSide2D::RightHandSide2D( const std::string& name ) :
  solver::TermBase<2,4,16,0>(name)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("background",m_background).link_to(&m_background)
      .mark_basic()
      .description("background flow");
  options().add("background_gradient",m_background_gradient).link_to(&m_background_gradient)
      .mark_basic()
      .description("background flow gradient");
  options().add("bdry_background",m_bdry_background).link_to(&m_bdry_background)
      .mark_basic()
      .description("boundary background flow");

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
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += L * solution()->array()[nodes[n]][eq];
    }
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[NEQS+eq] += L * m_background->array()[nodes[n]][eq];
    }
    if (ENABLE_SOURCE)
    {
      for (Uint eq=0; eq<NEQS; ++eq)
      {
        for (Uint d=0; d<NDIM; ++d)
        {
          vars[2*NEQS+NDIM*eq+d] += L * m_background_gradient->array()[nodes[n]][NDIM*eq+d];
        }
      }
    }
  }
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
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq]      += L * bdry_solution()->array()[nodes[n]][eq];
    }
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[NEQS+eq] += L * m_bdry_background->array()[nodes[n]][eq];
    }

  }
  //vars[NEQS+0] = m_rho0;
  //vars[NEQS+1] = m_U0[XX];
  //vars[NEQS+2] = m_U0[YY];
  //vars[NEQS+3] = m_p0;

  for (Uint eq=0; eq<NEQS; ++eq)
  {
    for (Uint d=0; d<NDIM; ++d)
    {
      vars[2*NEQS+NDIM*eq+d] = 0.;
    }
  }
}


void RightHandSide2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
}

void RightHandSide2D::compute_phys_data( const ColVector_NDIM& coords,
                                const RowVector_NVAR& vars,
                                const RowVector_NGRAD& gradvars,
                                const Matrix_NDIMxNGRAD& gradvars_grad,
                                DATA& phys_data )
{
  phys_data.rho0   = vars[NEQS+0];
  phys_data.U0[XX] = vars[NEQS+1];
  phys_data.U0[YY] = vars[NEQS+2];
  phys_data.p0     = vars[NEQS+3];
  phys_data.c0=std::sqrt(phys_data.gamma*phys_data.p0/phys_data.rho0);

  cf3_always_assert_desc(common::to_str(phys_data.rho0),std::abs(phys_data.rho0-1.) < 1e-7);
  cf3_always_assert_desc(common::to_str(phys_data.p0),std::abs(phys_data.p0-1.) < 1e-7);

  if( ENABLE_SOURCE )
  {
    for(Uint d=0; d<NDIM; ++d)
    {
      phys_data.grad_rho0[d] = vars[2*NEQS+d+NDIM*0];
      phys_data.grad_u0[d]   = vars[2*NEQS+d+NDIM*1];
      phys_data.grad_v0[d]   = vars[2*NEQS+d+NDIM*2];
      phys_data.grad_p0[d]   = vars[2*NEQS+d+NDIM*3];
    }
  }
  phys_data.compute_from_conservative(vars.head<4>());
}
void RightHandSide2D::compute_source( const DATA& phys_data, RowVector_NEQS& source )
{
  ///      [                                0                                   ]
  ///  H = [    (rho' u0 + rho0 u') du0/dx  +  (rho' v0 + rho0 v') du0/dy       ]
  ///      [    (rho' u0 + rho0 u') dv0/dx  +  (rho' v0 + rho0 v') dv0/dy       ]
  ///      [ (gamma-1) p' (du0/dx + dv0/dy) - (gamma-1) (u' dp0/dx + v' dp0/dy) ]

  const Real&       rho0    = phys_data.rho0;
  const Real&       u0      = phys_data.U0[XX];
  const Real&       v0      = phys_data.U0[YY];
  const Real&       p0      = phys_data.p0;
  const Real&       rho     = phys_data.rho;  // rho'
  const Real        rho0u   = phys_data.rho0 * phys_data.U[XX];  // rho0 u'
  const Real        rho0v   = phys_data.rho0 * phys_data.U[YY];  // rho0 v'
  const Real&       p       = phys_data.p;  // p'
  const Real        gm1     = m_gamma-1.;
  const Real&       u       = phys_data.U[XX];
  const Real&       v       = phys_data.U[YY];
  const ColVector_NDIM& grad_rho0 = phys_data.grad_rho0;
  const ColVector_NDIM& grad_u0   = phys_data.grad_u0;
  const ColVector_NDIM& grad_v0   = phys_data.grad_v0;
  const ColVector_NDIM& grad_p0   = phys_data.grad_p0;

  source[0] = - (                                    0                                    );
  source[1] = - (       (rho*u0 + rho0u)*grad_u0[XX] + (rho*v0 + rho0v)*grad_u0[YY]       );
  source[2] = - (       (rho*u0 + rho0u)*grad_v0[XX] + (rho*v0 + rho0v)*grad_v0[YY]       );
  source[3] = - ( gm1*p*(grad_u0[XX] + grad_v0[YY]) - gm1*(u*grad_p0[XX] + v*grad_p0[YY]) );
  // minus signs in front because this definition assumed H is on left hand side of equation, while this will
  // be added to right hand side
}


////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3
