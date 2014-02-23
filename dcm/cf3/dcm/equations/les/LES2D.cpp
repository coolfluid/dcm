// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/les/LES2D.hpp"
#include "cf3/solver/Term.hpp"

#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/ShapeFunction.hpp"
#include "cf3/mesh/ElementType.hpp"
#include "cf3/mesh/Connectivity.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < LES2D, solver::PDE, LibLES > LES2D_Builder;

////////////////////////////////////////////////////////////////////////////////

LES2D::LES2D ( const std::string& name  ) :
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
  options().add("kappa",2.601e-2)
      .description("Thermal conduction")
      .mark_basic();
  options().add("mu",1.806e-5)
      .description("Dynamic viscosity")
      .mark_basic();
  options().add("riemann_solver",std::string("Rusanov"))
      .mark_basic()
      .description("Riemann Solver");
  options().add("sfs_model",std::string("WALE"))
      .mark_basic()
      .description("Sub-filter-scale model");


  m_nb_dim = 2;
  m_nb_eqs = 4;

  add_time();
}

////////////////////////////////////////////////////////////////////////////////

LES2D::~LES2D()
{
}

////////////////////////////////////////////////////////////////////////////////

std::string LES2D::solution_variables() const
{
  return "rho[scalar], rhoU[vector], rhoE[scalar]";
}

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3
