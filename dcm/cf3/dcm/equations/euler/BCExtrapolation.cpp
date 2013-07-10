// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/euler/BCExtrapolation.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCExtrapolation1D,solver::BC,LibEuler> BCExtrapolation1D_builder;
common::ComponentBuilder<BCExtrapolation2D,solver::BC,LibEuler> BCExtrapolation2D_builder;
common::ComponentBuilder<BCExtrapolation3D,solver::BC,LibEuler> BCExtrapolation3D_builder;

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
