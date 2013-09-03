// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/euler/RiemannSolvers.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<Roe1D,common::Component,LibEuler> Roe1D_builder;
common::ComponentBuilder<Roe2D,common::Component,LibEuler> Roe2D_builder;

common::ComponentBuilder<HLLE1D,common::Component,LibEuler> HLLE1D_builder;
common::ComponentBuilder<HLLE2D,common::Component,LibEuler> HLLE2D_builder;

common::ComponentBuilder<Rusanov1D,common::Component,LibEuler> Rusanov1D_builder;
common::ComponentBuilder<Rusanov2D,common::Component,LibEuler> Rusanov2D_builder;

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
