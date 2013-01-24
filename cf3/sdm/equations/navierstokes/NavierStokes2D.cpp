// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/sdm/equations/navierstokes/NavierStokes2D.hpp"
#include "cf3/solver/Term.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < NavierStokes2D, solver::PDE, LibNavierStokes > NavierStokes2D_Builder;

////////////////////////////////////////////////////////////////////////////////

NavierStokes2D::NavierStokes2D ( const std::string& name  ) :
  PDE ( name )
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
  options().add("k",2.601e-2)
      .description("Heat conduction")
      .mark_basic();
  options().add("mu",1.806e-5)
      .description("Dynamic viscosity")
      .mark_basic();
  options().add("riemann_solver",std::string("Roe"))
      .mark_basic()
      .description("Riemann Solver");

  m_nb_dim = 2;
  m_nb_eqs = 4;

  add_time();
  add_term("terms","cf3.sdm.equations.navierstokes.NSTerms2D");
}

////////////////////////////////////////////////////////////////////////////////

NavierStokes2D::~NavierStokes2D()
{
}

////////////////////////////////////////////////////////////////////////////////

std::string NavierStokes2D::solution_variables() const
{
  return "rho[scalar], U[vector], rhoE[scalar]";
}

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3
