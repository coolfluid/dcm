// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/sdm/core/ShapeFunction.hpp"

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

ShapeFunction::ShapeFunction(const std::string& name) : mesh::ShapeFunction(name)
{
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3
