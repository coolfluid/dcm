// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/RegistLibrary.hpp"

#include "cf3/dcm/solver/erk/LibERK.hpp"

namespace cf3 {
namespace dcm {
namespace solver {
namespace erk {
  
using namespace common;

cf3::common::RegistLibrary<LibERK> libExplicit;

////////////////////////////////////////////////////////////////////////////////

void LibERK::initiate()
{
  if(m_is_initiated)
    return;

  Handle<Component> lib = Core::instance().libraries().get_child("cf3.dcm.solver.erk");
  cf3_assert(lib);

  m_is_initiated = true;
}

////////////////////////////////////////////////////////////////////////////////

} // erk
} // solver
} // dcm
} // cf3
