// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/dcm/equations/navierstokes/ComputePrimitiveVariables.hpp"
#include "cf3/mesh/Field.hpp"

using namespace cf3::common;
using namespace cf3::mesh;

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ComputePrimitiveVariables, common::Action, LibNavierStokes > ComputePrimitiveVariables_Builder;

////////////////////////////////////////////////////////////////////////////////

ComputePrimitiveVariables::ComputePrimitiveVariables ( const std::string& name  ) :
  common::Action ( name )
{
  // properties
  properties()["brief"] = std::string("NavierStokes 2D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 2D NavierStokes physics right-hand-side");

  options().add("gamma",1.4)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("R",287.05)
      .mark_basic()
      .description("Gas constant");
  options().add("solution",m_solution)
      .description("Solution field")
      .mark_basic()
      .link_to(&m_solution);
  options().add("pressure",m_pressure)
      .description("Pressure field")
      .mark_basic()
      .link_to(&m_pressure);
  options().add("velocity",m_velocity)
      .description("Velocity field")
      .mark_basic()
      .link_to(&m_velocity);
  options().add("temperature",m_temperature)
      .description("Temperature field")
      .mark_basic()
      .link_to(&m_temperature);
}

////////////////////////////////////////////////////////////////////////////////

ComputePrimitiveVariables::~ComputePrimitiveVariables()
{
}

////////////////////////////////////////////////////////////////////////////////

void ComputePrimitiveVariables::execute()
{
  Dictionary& dictionary = m_solution->dict();

  m_pressure = Handle<mesh::Field>( dictionary.get_child("pressure") );
  if( is_null(m_pressure) )
  {
    m_pressure = dictionary.create_field("pressure","p").handle<Field>();
    options().set("pressure",m_pressure);
  }

  m_velocity = Handle<mesh::Field>( dictionary.get_child("velocity") );
  if( is_null(m_velocity) )
  {
    m_velocity = dictionary.create_field("velocity","U[vector]").handle<Field>();
    options().set("velocity",m_velocity);
  }

  m_temperature = Handle<mesh::Field>( dictionary.get_child("temperature") );
  if( is_null(m_temperature) )
  {
    m_temperature = dictionary.create_field("temperature","T").handle<Field>();
    options().set("temperature",m_temperature);
  }

  const Uint nb_nodes = m_solution->size();
  const Uint dim = m_solution->dict().options().value<Uint>("dimension");
  const Real g = options().value<Real>("gamma");
  const Real R = options().value<Real>("R");
  Field::ArrayT& cons = m_solution->array();
  Field::ArrayT& p = m_pressure->array();
  Field::ArrayT& U = m_velocity->array();
  Field::ArrayT& T = m_temperature->array();

  for (Uint n=0; n<nb_nodes; ++n)
  {
    const Real rho = cons[n][0];
    const Real inv_rho = 1./rho;
    Real U2(0.);
    for (Uint d=0; d<dim; ++d)
    {
      U[n][d] = cons[n][d+1] * inv_rho;
      U2 += U[n][d]*U[n][d];
    }
    const Real rhoE = cons[n][dim+1];
    p[n][0] = (g-1.)*( rhoE - 0.5*rho*U2 );
    T[n][0] = p[n][0]/(rho*R);
  }
}

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3
