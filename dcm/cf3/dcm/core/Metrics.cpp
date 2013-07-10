// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/core/Metrics.hpp"

namespace cf3 {
namespace dcm {
namespace core {

///////////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < Metrics<1u>, common::Component, LibCore > Metrics1d_Builder;
common::ComponentBuilder < Metrics<2u>, common::Component, LibCore > Metrics2d_Builder;
common::ComponentBuilder < Metrics<3u>, common::Component, LibCore > Metrics3d_Builder;

/////////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3
