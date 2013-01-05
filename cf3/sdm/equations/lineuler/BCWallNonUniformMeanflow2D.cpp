// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/sdm/equations/lineuler/BCWallNonUniformMeanflow2D.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCWallNonUniformMeanflow2D,BC,LibLinEuler>BCWallNonUniformMeanflow2D_builder;

/////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
