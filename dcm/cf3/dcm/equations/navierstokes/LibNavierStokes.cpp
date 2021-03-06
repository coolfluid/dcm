// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/RegistLibrary.hpp"

#include "cf3/dcm/equations/navierstokes/LibNavierStokes.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

  using namespace common;

cf3::common::RegistLibrary<LibNavierStokes> LibNavierStokes;


} // navierstokes
} // equations
} // dcm
} // cf3

