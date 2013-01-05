// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/sdm/equations/euler/Rusanov.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<Rusanov1D,solver::RiemannSolver<Terms1D>,LibEuler> Rusanov1D_builder;
common::ComponentBuilder<Rusanov2D,solver::RiemannSolver<Terms2D>,LibEuler> Rusanov2D_builder;

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3
