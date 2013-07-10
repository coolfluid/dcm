// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/navierstokes/BCWall.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCWall1D,solver::BC,LibNavierStokes> BCWall1D_builder;
common::ComponentBuilder<BCWall2D,solver::BC,LibNavierStokes> BCWall2D_builder;

////////////////////////////////////////////////////////////////////////////////

void BCWall1D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
{
  boundary_solution[0] =  inner_solution[0];
  boundary_solution[1] = -inner_solution[1];
  boundary_solution[2] =  inner_solution[2];
}

////////////////////////////////////////////////////////////////////////////////

BCWall2D::BCWall2D(const std::string& name) :
  dcm::core::BC<2,4>(name),
  m_wall_velocity(0.),
  m_R(0.),
  m_gamma(0.)
{
  options().add("wall_velocity",m_wall_velocity)
      .description("The velocity of the wall")
      .link_to(&m_wall_velocity)
      .mark_basic();

  options().add("gamma",m_gamma).link_to(&m_gamma);
  options().add("R",m_R).link_to(&m_R);
}
void BCWall2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
{
  boundary_solution[0]=  inner_solution[0];
  boundary_solution[1]= -inner_solution[1] - 2.*inner_solution[0]*face_normal[YY]*m_wall_velocity;
  boundary_solution[2]= -inner_solution[2] + 2.*inner_solution[0]*face_normal[XX]*m_wall_velocity;
  boundary_solution[3]=  inner_solution[3];

  boundary_solution_gradient.col(0) = inner_solution_gradient.col(0);
  boundary_solution_gradient.col(1) = inner_solution_gradient.col(1);
  boundary_solution_gradient.col(2) = inner_solution_gradient.col(2);

  // Set zero temperature gradient normal to the wall
  cf3_assert(m_gamma!=0);
  cf3_assert(m_R!=0);
  const Real rho = inner_solution[0];
  const Real u = inner_solution[1]/rho;
  const Real v = inner_solution[2]/rho;
  const Real E = inner_solution[3]/rho;
  const Real U2 = u*u+v*v;
  const Real p  = (m_gamma-1.)*rho*(E - 0.5*U2);
  const Real T  = p/(rho*m_R);

  const ColVector_NDIM grad_rho  = inner_solution_gradient.col(0);
  const ColVector_NDIM grad_rhou = inner_solution_gradient.col(1);
  const ColVector_NDIM grad_rhov = inner_solution_gradient.col(2);
  const ColVector_NDIM grad_rhoE = inner_solution_gradient.col(3);

  const ColVector_NDIM grad_u = 1./rho * ( grad_rhou - u*grad_rho );
  const ColVector_NDIM grad_v = 1./rho * ( grad_rhov - v*grad_rho );
  const ColVector_NDIM grad_T = ( (grad_rhoE - u*grad_u - v*grad_v)*(m_gamma-1)*rho -p*grad_rho )/(m_R*rho*rho);

  const Real dTdn = grad_T.dot(face_normal);
  const ColVector_NDIM bdry_grad_T = grad_T - 2.*dTdn*face_normal;

  const Real bdry_u = boundary_solution[1]/rho;
  const Real bdry_v = boundary_solution[2]/rho;
  // transform to gradient of rhoE
  const ColVector_NDIM bdry_grad_rhoE =
      bdry_u*grad_u + bdry_v*grad_v + ( m_R*(T*grad_rho + rho*bdry_grad_T ) )/(m_gamma-1);
  boundary_solution_gradient.col(3) = bdry_grad_rhoE;
}

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3
