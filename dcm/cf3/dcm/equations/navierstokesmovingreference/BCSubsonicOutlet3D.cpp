// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "common/Builder.hpp"

#include "dcm/navierstokesmovingreference/BCSubsonicOutlet3D.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace navierstokesmovingreference {
	
//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCSubsonicOutlet3D,BC,LibNavierStokesMovingReference> bcsubsonicoutlet3d_builder;

/////////////////////////////////////////////////////////////////////////////

} // navierstokesmovingreference
} // dcm
} // cf3
