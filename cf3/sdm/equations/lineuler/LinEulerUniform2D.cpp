// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/sdm/equations/lineuler/LinEulerUniform2D.hpp"
#include "cf3/solver/Term.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::solver;

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < LinEulerUniform2D, solver::PDE, LibLinEuler > LinEulerUniform2D_Builder;

////////////////////////////////////////////////////////////////////////////////

LinEulerUniform2D::LinEulerUniform2D ( const std::string& name  ) :
  PDE ( name )
{
  // properties
  properties()["brief"] = std::string("Linearized Euler 2D with uniform mean flow");
  properties()["description"] = std::string("Component that can solve the 1D Euler physics right-hand-side");

  options().add("gamma",1.4).mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("rho0",1.0).mark_basic()
      .description("Constant mean density");
  options().add("U0",std::vector<Real>(2,0.)).mark_basic()
      .description("Constant mean velocity");
  options().add("p0",1.).mark_basic()
      .description("Constant mean pressure");

  m_nb_dim = 2;
  m_nb_eqs = 4;

  add_time();
  add_term("terms","cf3.sdm.equations.lineuler.TermsUniform2D");
}

////////////////////////////////////////////////////////////////////////////////

LinEulerUniform2D::~LinEulerUniform2D()
{
}

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
