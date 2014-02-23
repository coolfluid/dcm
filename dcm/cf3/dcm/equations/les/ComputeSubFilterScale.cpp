// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/dcm/equations/les/ComputeSubFilterScale.hpp"
#include "cf3/dcm/equations/les/EddyViscosityModel.hpp"
#include "cf3/dcm/core/ShapeFunction.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/ShapeFunction.hpp"
#include "cf3/dcm/core/Reconstructions.hpp"
#include "cf3/dcm/core/Metrics.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::dcm::core;

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ComputeSubFilterScale, common::Action, LibLES > ComputeSubFilterScale_Builder;

////////////////////////////////////////////////////////////////////////////////

ComputeSubFilterScale::ComputeSubFilterScale ( const std::string& name  ) :
  common::Action ( name )
{
  // properties
  properties()["brief"] = std::string("NavierStokes 2D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 2D NavierStokes physics right-hand-side");

  options().add("sfs_model", m_sfs_model)
    .description("Sub-filter-scale model")
    .mark_basic()
    .link_to(&m_sfs_model);
  options().add("density",m_density)
      .description("Density field")
      .mark_basic()
      .link_to(&m_density);
  options().add("velocity",m_velocity)
      .description("Velocity field")
      .mark_basic()
      .link_to(&m_velocity);
  options().add("sfs_kinetic_energy",m_k_sfs)
      .description("[out] subfilterscale_kinetic_energy")
      .mark_basic()
      .link_to(&m_k_sfs);
  options().add("sfs_viscosity",m_nuT)
      .description("[out] subfilterscale_kinetic_energy")
      .mark_basic()
      .link_to(&m_nuT);
  options().add("sfs_heat_conduction",m_kappaT)
      .description("[out] subfilterscale_heat_conduction")
      .mark_basic()
      .link_to(&m_kappaT);
}

////////////////////////////////////////////////////////////////////////////////

ComputeSubFilterScale::~ComputeSubFilterScale()
{
}

////////////////////////////////////////////////////////////////////////////////

void ComputeSubFilterScale::execute()
{
  Dictionary& dictionary = m_velocity->dict();

  if( is_null(m_velocity) ) throw SetupError(FromHere(), "velocity field not setup");
  if( is_null(m_density)  ) throw SetupError(FromHere(), "density field not setup");

  const Uint dim = dictionary.options().value<Uint>("dimension");

  Dictionary& dict = m_velocity->dict();

  Real rho;
  RealVector2 U;
  RealMatrix2 grad_U;
  Real filter_width;
  Real Jdet;

  int NDIM=2;

  boost_foreach( const Handle<Space>& space, dict.spaces() )
  {
    const dcm::core::ShapeFunction& sf = *space->shape_function().handle<dcm::core::ShapeFunction>();
        
    Handle< dcm::core::Metrics<2u> > metrics( space->get_child("metrics") );
    if (is_null(metrics))
    {
      metrics = space->create_component< dcm::core::Metrics<2u> >("metrics");
      metrics->setup_for_space(space);
    }

    Uint nb_cells = space->size();
    Uint nb_nodes_per_cell = sf.nb_nodes();
    for (Uint cell_idx=0; cell_idx<nb_cells; ++cell_idx)
    {
      ElementMetrics<2u>& element_metrics = *metrics->element(cell_idx);
      mesh::Connectivity::ConstRow nodes = space->connectivity()[cell_idx];

      for( Uint pt=0; pt<nb_nodes_per_cell; ++pt)
      {
        Jdet = element_metrics.sol_pt_Jdet(pt);
        filter_width = std::pow(Jdet/nb_nodes_per_cell, 1./static_cast<Real>(NDIM));

        rho   = m_density->array()[  nodes[pt] ][0];
        U[XX] = m_velocity->array()[ nodes[pt] ][XX];
        U[XX] = m_velocity->array()[ nodes[pt] ][YY];
        grad_U.setZero();
        for (Uint d=0; d<2; ++d)
        {
          boost_foreach( const Uint sol_pt, metrics->gradient_from_sol_pts_to_sol_pt(pt)[d].used_points() )
          {
            const Real n = nodes[sol_pt];
            const Real D = metrics->gradient_from_sol_pts_to_sol_pt(pt)[d].coeff(sol_pt);
            grad_U(d,XX) += D * m_velocity->array()[n][XX];
            grad_U(d,YY) += D * m_velocity->array()[n][YY];
          }
        }
        grad_U = element_metrics.sol_pt_Jinv(pt) * grad_U;

        m_sfs_model->compute( rho, U, grad_U, filter_width );

        m_k_sfs->array()[  nodes[pt] ][0] = m_sfs_model->sfs_kinetic_energy();
        m_nuT->array()[    nodes[pt] ][0] = m_sfs_model->sfs_viscosity();
        m_kappaT->array()[ nodes[pt] ][0] = m_sfs_model->sfs_heat_conduction();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3
