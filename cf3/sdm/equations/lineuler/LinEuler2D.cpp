// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "cf3/common/Builder.hpp"
#include "cf3/common/Log.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/Group.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/common/OptionArray.hpp"
#include "cf3/common/Signal.hpp"
#include "cf3/math/Defs.hpp"
#include "cf3/math/Consts.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Entities.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/sdm/BoundaryConditions.hpp"
#include "cf3/sdm/Physics.hpp"
#include "cf3/sdm/DomainDiscretization.hpp"
#include "cf3/sdm/Term.hpp"
#include "cf3/sdm/equations/lineuler/LinEuler2D.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::physics;
using namespace cf3::solver;


namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < LinEuler2D, sdm::Physics, LibLinEuler > LinEuler2D_Builder;

////////////////////////////////////////////////////////////////////////////////

LinEuler2D::LinEuler2D ( const std::string& name  ) :
  sdm::Physics ( name )
{
  // properties

  properties()["brief"] = std::string("Linearized Euler Equations with Non-Uniform mean flow");
  properties()["description"] = std::string("Long description not available");
  
  Uint dim = 2;
  m_nb_eqs = 4u;
  
  m_gamma = 1.4;
  options().add("gamma",m_gamma).link_to(&m_gamma).mark_basic();
}

////////////////////////////////////////////////////////////////////////////////

LinEuler2D::~LinEuler2D()
{
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::create_terms( const Handle<sdm::DomainDiscretization>& domain_discretization)
{
  m_convection = domain_discretization->create_term("convection","cf3.sdm.lineuler.ConvectionNonUniformMeanflow2D").handle<Term>();
  configure_terms();
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::create_fields( const Handle<mesh::Dictionary>& dict )
{
  m_coords    = dict->coordinates().handle<Field>();
  m_solution  = dict->create_field("solution" , "rho[s],rho0_U[vector],p[s]" ).handle<Field>();

  m_rho0 = dict->create_field("rho0", SCALAR    ).handle<Field>();
  m_U0   = dict->create_field("U0",   VECTOR_2D ).handle<Field>();
  m_p0   = dict->create_field("p0",   SCALAR    ).handle<Field>();

  m_grad_U0 = dict->create_field("grad_U0", TENSOR_2D ).handle<Field>();
  m_grad_p0 = dict->create_field("grad_p0", VECTOR_2D ).handle<Field>();

  configure_terms();
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::configure_terms()
{
  if (m_convection)
  {
    m_convection->options().set("gamma", options().value<Real>("gamma") );
    if (m_rho0)
      m_convection->options().set("rho0", m_rho0 );
    if (m_U0)
      m_convection->options().set("U0", m_U0 );
    if (m_p0)
      m_convection->options().set("p0", m_p0 );
    if (m_solution)
      m_convection->options().set("solution", m_solution);
  }  
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_sol_pt_phys_data( const mesh::SpaceElem& cell, const Uint sol_pt,
                                           PhysData& phys_data )
{
  const Uint n = cell.nodes()[sol_pt];

  phys_data.coord[XX] = m_coords->array()[n][XX];
  phys_data.coord[YY] = m_coords->array()[n][YY];

  phys_data.solution[0]    = m_solution->array()[n][0];
  phys_data.solution[1]    = m_solution->array()[n][1];
  phys_data.solution[2]    = m_solution->array()[n][2];
  phys_data.solution[3]    = m_solution->array()[n][3];

  phys_data.rho0      = m_rho0->array()[n][0];
  phys_data.U0[XX]    = m_U0->array()[n][XX];
  phys_data.U0[YY]    = m_U0->array()[n][YY];
  phys_data.p0        = m_p0->array()[n][0];

//  phys_data.grad_U0(XX,XX) = m_grad_U0->array()[n][0];
//  phys_data.grad_U0(YY,XX) = m_grad_U0->array()[n][1];
//  phys_data.grad_U0(XX,YY) = m_grad_U0->array()[n][2];
//  phys_data.grad_U0(YY,YY) = m_grad_U0->array()[n][3];
//  phys_data.grad_p0[XX]    = m_grad_p0->array()[n][XX];
//  phys_data.grad_p0[YY]    = m_grad_p0->array()[n][YY];

  // computations
  phys_data.gamma = m_gamma;
  phys_data.compute();
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_flx_pt_phys_data( const mesh::SpaceElem& cell, const Uint flx_pt,
                                           PhysData& phys_data )
{
  const ReconstructPoint& solution_to_flux = m_solution_to_flux[cell.comp][flx_pt];
  const ReconstructPoint& geometry_to_flux = m_geometry_to_flux[cell.comp][flx_pt];

  RealMatrix coords = cell.comp->support().geometry_space().get_coordinates(cell.idx);
  geometry_to_flux(coords , phys_data.coord);

  enum InterpolationVars { sol0, sol1, sol2, sol3, rho0, u0, v0, p0 };
  RealMatrix in_solution_pts(cell.shape_function().nb_nodes(),8u);
  RealVector from_solution_pts(8u);

  // Only loop over solution points that will contribute to this flux point for efficiency
  boost_foreach( Uint sol_pt, solution_to_flux.used_points() )
  {
    const Uint n = cell.nodes()[sol_pt];
    in_solution_pts(sol_pt,sol0)  = m_solution->array()[n][0];
    in_solution_pts(sol_pt,sol1)  = m_solution->array()[n][1];
    in_solution_pts(sol_pt,sol2)  = m_solution->array()[n][2];
    in_solution_pts(sol_pt,sol3)  = m_solution->array()[n][3];

    in_solution_pts(sol_pt,rho0)  = m_rho0->array()[n][0];
    in_solution_pts(sol_pt,u0)    = m_U0->array()[n][XX];
    in_solution_pts(sol_pt,v0)    = m_U0->array()[n][YY];
    in_solution_pts(sol_pt,p0)    = m_p0->array()[n][0];
  }
  solution_to_flux(in_solution_pts, from_solution_pts);

  phys_data.solution[0]  = from_solution_pts[sol0];
  phys_data.solution[1]  = from_solution_pts[sol1];
  phys_data.solution[2]  = from_solution_pts[sol2];
  phys_data.solution[3]  = from_solution_pts[sol3];
  phys_data.rho0    = from_solution_pts[rho0];
  phys_data.U0[XX]  = from_solution_pts[u0];
  phys_data.U0[YY]  = from_solution_pts[v0];
  phys_data.p0      = from_solution_pts[p0];

  // computations
  phys_data.gamma = m_gamma;
  phys_data.compute();

}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_convective_flux(const PhysData& p, const RealVectorNDIM &normal,
                                         RealVectorNEQS &flux)
{
  const Real u0n = p.U0.dot(normal);
  const Real un  = p.U.dot(normal);

  flux[0] = u0n*p.rho        + p.rho0*un;
  flux[1] = u0n*p.rho0_U[XX] + p.p*normal[XX];
  flux[2] = u0n*p.rho0_U[YY] + p.p*normal[YY];
  flux[3] = u0n*p.p          + p.rho0*un*p.c0*p.c0;
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_convective_wave_speed( const PhysData& p, const RealVectorNDIM& normal,
                                               Real& wave_speed)
{
  const Real u0n = p.U0.dot(normal);
  wave_speed = std::abs(u0n)+p.c0;
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_convective_eigenvalues( const PhysData& p, const RealVectorNDIM& normal,
                                                 RealVectorNEQS& eigenvalues)
{
  const Real u0n = p.U0.dot(normal);

  eigenvalues <<
      u0n,
      u0n,
      u0n + p.c0,
      u0n - p.c0;
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_convective_eigensystem( const PhysData& p, const RealVectorNDIM& normal,
                                                 RealMatrixNEQSxNEQS& right_eigenvectors,
                                                 RealVectorNEQS& eigenvalues,
                                                 RealMatrixNEQSxNEQS& left_eigenvectors)
{
  const Real u0n = p.U0.dot(normal);

  const Real  inv_c  = 1./p.c0;
  const Real  inv_c2 = inv_c*inv_c;
  const Real& nx = normal[XX];
  const Real& ny = normal[YY];

  right_eigenvectors <<
        1.,      0.,      0.5*inv_c,    0.5*inv_c,
        0.,      ny,      0.5*nx,      -0.5*nx,
        0.,     -nx,      0.5*ny,      -0.5*ny,
        0.,      0.,      0.5*p.c0,     0.5*p.c0;

  left_eigenvectors <<
        1.,      0.,      0.,          -inv_c2,
        0.,      ny,     -nx,           0,
        0.,      nx,      ny,           inv_c,
        0.,     -nx,     -ny,           inv_c;

  eigenvalues <<
      u0n,
      u0n,
      u0n + p.c0,
      u0n - p.c0;
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_absolute_flux_jacobian( const PhysData& p, const RealVectorNDIM& normal,
                                                 RealMatrixNEQSxNEQS& absolute_flux_jacobian)
{
  const Real u0n = p.U0.dot(normal);
  const Real inv_2c  = 0.5/p.c0;
  const Real inv_2c2 = 0.5/(p.c0*p.c0);

  const Real& nx = normal[XX];
  const Real& ny = normal[YY];
  const Real nx2 = nx*nx;
  const Real ny2 = ny*ny;
  const Real absu0n = std::abs(u0n);
  const Real cpu = std::abs(p.c0+u0n);
  const Real cmu = std::abs(p.c0-u0n);
  const Real plus  = cmu + cpu;
  const Real minus = cpu - cmu;
  const Real pm2u  = plus - 2*absu0n;

  absolute_flux_jacobian <<
    absu0n,     (nx*minus)*inv_2c,              (ny*minus)*inv_2c,              pm2u*inv_2c2,
    0,          (2*ny2*absu0n + nx2*plus)*0.5,  (nx*ny*pm2u)*0.5,               (nx*minus)*inv_2c,
    0,          (nx*ny*pm2u)*0.5,               (2*nx2*absu0n + ny2*plus)*0.5,  (ny*minus)*inv_2c,
    0,          (p.c0*nx*minus)*0.5,            (p.c0*ny*minus)*0.5,            plus*0.5;
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_interior_fluxes(const PhysData& p, const RealVectorNDIM& normal,
                                          RealVectorNEQS& flux, Real &wave_speed)
{
  compute_convective_flux(p,normal,flux);
  compute_convective_wave_speed(p,normal,wave_speed);
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::compute_face_fluxes(const PhysData& left, const PhysData& right, const RealVectorNDIM& normal,
                                     RealVectorNEQS& flux, Real& wave_speed)
{
  compute_convective_flux(left, normal,flux_left);
  compute_convective_flux(right,normal,flux_right);

  // No Roe-average is needed as the eigen-system is only dependant of the mean flow
  compute_absolute_flux_jacobian(left,normal,A);

  flux.noalias()  = 0.5*(flux_left+flux_right);
  flux.noalias() -= 0.5*A*(right.solution-left.solution);

  compute_convective_wave_speed(left,normal,wave_speed);
}

////////////////////////////////////////////////////////////////////////////////

LinEuler2D::RealMatrixNEQSxNEQS LinEuler2D::A;
LinEuler2D::RealMatrixNEQSxNEQS LinEuler2D::L;
LinEuler2D::RealMatrixNEQSxNEQS LinEuler2D::R;
LinEuler2D::RealVectorNEQS LinEuler2D::D;
LinEuler2D::RealVectorNEQS LinEuler2D::flux_left;
LinEuler2D::RealVectorNEQS LinEuler2D::flux_right;

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
