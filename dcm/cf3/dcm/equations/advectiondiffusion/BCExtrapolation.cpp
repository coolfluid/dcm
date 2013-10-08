// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/advectiondiffusion/BCExtrapolation.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCExtrapolation<1>,solver::BC,LibAdvectionDiffusion> BCExtrapolation1D_builder;
common::ComponentBuilder<BCExtrapolation<2>,solver::BC,LibAdvectionDiffusion> BCExtrapolation2D_builder;
common::ComponentBuilder<BCExtrapolation<3>,solver::BC,LibAdvectionDiffusion> BCExtrapolation3D_builder;

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3
