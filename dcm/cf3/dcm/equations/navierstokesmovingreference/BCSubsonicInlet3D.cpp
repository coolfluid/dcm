// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "common/Builder.hpp"

#include "dcm/navierstokesmovingreference/BCSubsonicInlet3D.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace navierstokesmovingreference {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCSubsonicInletUT3D,BC,LibNavierStokesMovingReference> bcsubsonicinletUT3d_builder;
common::ComponentBuilder<BCSubsonicInletTtPtAlpha3D,BC,LibNavierStokesMovingReference> bcsubsonicinletTtPtAlpha3d_builder;

/////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // dcm
} // cf3
