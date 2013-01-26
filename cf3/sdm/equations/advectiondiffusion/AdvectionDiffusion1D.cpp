// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"

#include "cf3/sdm/equations/advectiondiffusion/AdvectionDiffusion1D.hpp"
#include "cf3/solver/Term.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace sdm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < AdvectionDiffusion1D, solver::PDE, LibAdvectionDiffusion > AdvectionDiffusion1D_Builder;

////////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion1D::AdvectionDiffusion1D ( const std::string& name  ) :
  PDE ( name )
{
  // properties
  properties()["brief"] = std::string("advectiondiffusion 1D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 1D advectiondiffusion physics right-hand-side");

  options().add("a",1.)
      .mark_basic()
      .description("Advection speed");
  options().add("mu",0.)
      .mark_basic()
      .description("Diffusion coefficient");

  m_nb_dim = 1;
  m_nb_eqs = 1;

  add_time();
  add_term("terms","cf3.sdm.equations.advectiondiffusion.Terms1D");
}

////////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion1D::~AdvectionDiffusion1D()
{
}

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // sdm
} // cf3
