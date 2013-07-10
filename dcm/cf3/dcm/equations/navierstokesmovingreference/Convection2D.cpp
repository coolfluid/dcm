// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "common/Builder.hpp"

#include "dcm/navierstokesmovingreference/Convection2D.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace navierstokesmovingreference {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<Convection2D,Term,LibNavierStokesMovingReference> convection2d_builder;

/////////////////////////////////////////////////////////////////////////////

} // navierstokesmovingreference
} // dcm
} // cf3
