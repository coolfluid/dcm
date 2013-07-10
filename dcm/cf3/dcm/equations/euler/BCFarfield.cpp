// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/euler/BCFarfield.hpp"

using namespace cf3::common;

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCFarfield1D,solver::BC,LibEuler> BCFarfield1D_builder;
common::ComponentBuilder<BCFarfield2D,solver::BC,LibEuler> BCFarfield2D_builder;
common::ComponentBuilder<BCFarfield3D,solver::BC,LibEuler> BCFarfield3D_builder;

////////////////////////////////////////////////////////////////////////////////

BCFarfield1D::BCFarfield1D(const std::string& name)
   : dcm::core::BC<1,3>(name),
     m_rho(1.),
     m_p(1.),
     m_u(0.)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("rho",m_rho)
      .link_to(&m_rho)
      .mark_basic()
      .description("Farfield density");
  options().add("u",m_u)
      .link_to(&m_u)
      .mark_basic()
      .description("Farfield velocity");
  options().add("p",m_p)
      .link_to(&m_p)
      .mark_basic()
      .description("Farfield pressure");
  options()["rho"].attach_trigger( boost::bind( &BCFarfield1D::compute_state, this ) );
  options()["u"  ].attach_trigger( boost::bind( &BCFarfield1D::compute_state, this ) );
  options()["p"  ].attach_trigger( boost::bind( &BCFarfield1D::compute_state, this ) );
}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield1D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                              const ColVector_NDIM& coords,
                                              const ColVector_NDIM& face_normal, 
                                              RowVector_NEQS& boundary_solution )
{
  boundary_solution = state;
}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield1D::compute_state()
{
  state[0] = m_rho;
  state[1] = m_rho*m_u;
  state[2] = m_p/(m_gamma-1.) + 0.5*m_rho*m_u*m_u;
}

////////////////////////////////////////////////////////////////////////////////

BCFarfield2D::BCFarfield2D(const std::string& name) 
  : dcm::core::BC<2,4>(name),
    m_rho(1.),
    m_p(1.),
    m_U(2,0.)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("rho",m_rho)
      .link_to(&m_rho)
      .mark_basic()
      .description("Farfield density");
  options().add("U",m_U)
      .link_to(&m_U)
      .mark_basic()
      .description("Farfield velocity");
  options().add("p",m_p)
      .link_to(&m_p)
      .mark_basic()
      .description("Farfield pressure");
  options()["rho"].attach_trigger( boost::bind( &BCFarfield2D::compute_state, this ) );
  options()["U"  ].attach_trigger( boost::bind( &BCFarfield2D::compute_state, this ) );
  options()["p"  ].attach_trigger( boost::bind( &BCFarfield2D::compute_state, this ) );

}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield2D::compute_state()
{
  if (m_U.size()!=NDIM) throw SetupError(FromHere(), "U ["+to_str(m_U)+"] has wrong dimensions");
  state[0] = m_rho;
  state[1] = m_rho*m_U[XX];
  state[2] = m_rho*m_U[YY];
  state[3] = m_p/(m_gamma-1.) + 0.5*m_rho*(m_U[XX]*m_U[XX]+m_U[YY]*m_U[YY]);
}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield2D::compute_boundary_solution( const RowVector_NEQS& inner_solution, 
                                              const ColVector_NDIM& coords,
                                              const ColVector_NDIM& face_normal, 
                                              RowVector_NEQS& boundary_solution )
{
  boundary_solution = state;
}

////////////////////////////////////////////////////////////////////////////////

BCFarfield3D::BCFarfield3D(const std::string& name) 
  : dcm::core::BC<3,5>(name),
    m_rho(1.),
    m_p(1.),
    m_U(3,0.)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("rho",m_rho)
      .link_to(&m_rho)
      .mark_basic()
      .description("Farfield density");
  options().add("U",m_U)
      .link_to(&m_U)
      .mark_basic()
      .description("Farfield velocity");
  options().add("p",m_p)
      .link_to(&m_p)
      .mark_basic()
      .description("Farfield pressure");
  options()["rho"].attach_trigger( boost::bind( &BCFarfield3D::compute_state, this ) );
  options()["U"  ].attach_trigger( boost::bind( &BCFarfield3D::compute_state, this ) );
  options()["p"  ].attach_trigger( boost::bind( &BCFarfield3D::compute_state, this ) );
}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield3D::compute_state()
{
  if (m_U.size()!=NDIM) throw SetupError(FromHere(), "U has wrong dimensions");
  state[0] = m_rho;
  state[1] = m_rho*m_U[XX];
  state[2] = m_rho*m_U[YY];
  state[3] = m_rho*m_U[ZZ];
  state[4] = m_p/(m_gamma-1.) + 0.5*m_rho*(  m_U[XX]*m_U[XX]
                                            +m_U[YY]*m_U[YY]
                                            +m_U[ZZ]*m_U[ZZ] );
}

////////////////////////////////////////////////////////////////////////////////

void BCFarfield3D::compute_boundary_solution( const RowVector_NEQS& inner_solution, 
                                              const ColVector_NDIM& coords,
                                              const ColVector_NDIM& face_normal, 
                                              RowVector_NEQS& boundary_solution )
{
  boundary_solution = state;
}

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
