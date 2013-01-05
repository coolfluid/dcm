// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/RegistLibrary.hpp"

#include "cf3/sdm/solver/erkls/LibERKLS.hpp"

namespace cf3 {
namespace sdm {
namespace solver {
namespace erkls {
  
using namespace common;

cf3::common::RegistLibrary<LibERKLS> libERKLS;

////////////////////////////////////////////////////////////////////////////////

void LibERKLS::initiate()
{
  if(m_is_initiated)
    return;

  Handle<Component> lib = Core::instance().libraries().get_child("cf3.sdm.solver.erkls");
  cf3_assert(lib);

  m_is_initiated = true;
}

////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // sdm
} // cf3
