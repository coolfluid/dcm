// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/dcm/equations/lineuler/LinEuler2D.hpp"
#include "cf3/solver/Term.hpp"
#include "cf3/mesh/Field.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < LinEuler2D, solver::PDE, LibLinEuler > LinEuler2D_Builder;

////////////////////////////////////////////////////////////////////////////////

LinEuler2D::LinEuler2D ( const std::string& name  ) :
  PDE ( name )
{
  // properties
  properties()["brief"] = std::string("Linearized Euler 2D with uniform mean flow");
  properties()["description"] = std::string("Component that can solve the 1D Euler physics right-hand-side");

  options().add("gamma",1.).mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("background",m_background).link_to(&m_background).mark_basic()
      .description("Background flow");
  options().add("background_gradient",m_background_gradient).link_to(&m_background_gradient).mark_basic()
      .description("Gradient of background flow");
  options().add("bdry_background",m_bdry_background).link_to(&m_bdry_background).mark_basic()
      .description("Background flow in boundary");

  m_nb_dim = 2;
  m_nb_eqs = 4;

  add_time();
}

////////////////////////////////////////////////////////////////////////////////

LinEuler2D::~LinEuler2D()
{
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::create_fields()
{

  PDE::create_fields();
  CFinfo << "Creating LinEuler fields" << CFendl;
  if ( is_null(m_background) || ( &m_background->dict() != m_fields.get() ) )
  {
    if ( Handle<Component> found = m_fields->get_child("background") )
    {
      m_background = found->handle<Field>();
    }
    else
    {
      m_background = m_fields->create_field("background","rho0,U0[vec],p0").handle<Field>();
    }
  }
  options().set("background",m_background);

  if ( is_null(m_background_gradient) || ( &m_background_gradient->dict() != m_fields.get() ) )
  {
    if ( Handle<Component> found = m_fields->get_child("background_gradient") )
    {
      m_background_gradient = found->handle<Field>();
    }
    else
    {
      m_background_gradient = m_fields->create_field("background_gradient","grad_rho0[v],grad_u0[v],grad_v0[v],grad_p0[v]").handle<Field>();
      *m_background_gradient = 0.;
    }
  }
  options().set("background_gradient",m_background_gradient);
}

////////////////////////////////////////////////////////////////////////////////

void LinEuler2D::create_bdry_fields()
{

  PDE::create_bdry_fields();
  CFinfo << "Creating LinEuler bdry fields" << CFendl;

  if ( is_null(m_bdry_background) || ( &m_bdry_background->dict() != m_bdry_fields.get() ) )
  {
    if ( Handle<Component> found = m_bdry_fields->get_child("bdry_background") )
    {
      m_bdry_background = found->handle<Field>();
    }
    else
    {
      m_bdry_background = m_bdry_fields->create_field("bdry_background","rho0,U0[vec],p0").handle<Field>();
    }
  }
  options().set("bdry_background",m_bdry_background);

}

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3
