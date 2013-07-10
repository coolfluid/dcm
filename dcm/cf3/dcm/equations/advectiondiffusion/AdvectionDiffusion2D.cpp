// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"

#include "cf3/dcm/equations/advectiondiffusion/AdvectionDiffusion2D.hpp"
#include "cf3/solver/Term.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < AdvectionDiffusion2D, solver::PDE, LibAdvectionDiffusion > AdvectionDiffusion2D_Builder;

////////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion2D::AdvectionDiffusion2D ( const std::string& name  ) :
  PDE ( name )
{
  // properties
  properties()["brief"] = std::string("advectiondiffusion 2D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 2D advectiondiffusion physics right-hand-side");

  options().add("a",std::vector<Real>(2,0.))
      .mark_basic()
      .description("Advection speed");
  options().add("mu",0.)
      .mark_basic()
      .description("Diffusion coefficient");

  m_nb_dim = 2;
  m_nb_eqs = 1;

  add_time();
}

////////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion2D::~AdvectionDiffusion2D()
{
}

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3
