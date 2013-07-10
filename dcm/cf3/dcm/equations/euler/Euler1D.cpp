// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"

#include "cf3/dcm/equations/euler/Euler1D.hpp"
#include "cf3/solver/Term.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < Euler1D, solver::PDE, LibEuler > Euler1D_Builder;

////////////////////////////////////////////////////////////////////////////////

Euler1D::Euler1D ( const std::string& name  ) :
  PDE ( name )
{
  // properties
  properties()["brief"] = std::string("Euler 1D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 1D Euler physics right-hand-side");

  options().add("gamma",1.4)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("R",287.05)
      .mark_basic()
      .description("Gas constant");
  options().add("riemann_solver",std::string("Roe"))
      .mark_basic()
      .description("Riemann Solver");

  m_nb_dim = 1;
  m_nb_eqs = 3;

  add_time();
}

////////////////////////////////////////////////////////////////////////////////

Euler1D::~Euler1D()
{
}

////////////////////////////////////////////////////////////////////////////////

std::string Euler1D::solution_variables() const
{
  return "rho[scalar], rhoU[vector], rhoE[scalar]";
}

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
