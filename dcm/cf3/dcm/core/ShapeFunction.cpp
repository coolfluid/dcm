// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/dcm/core/ShapeFunction.hpp"

namespace cf3 {
namespace dcm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

ShapeFunction::ShapeFunction(const std::string& name) : mesh::ShapeFunction(name)
{
}

////////////////////////////////////////////////////////////////////////////////

const RealMatrix& ShapeFunction::mononomial_coefficients() const
{
  throw common::NotImplemented( FromHere(), "mononomial_coefficients not implemented");
  static RealMatrix dummy;
  return dummy;
}

////////////////////////////////////////////////////////////////////////////////

const RealMatrix& ShapeFunction::mononomial_exponents() const
{
  throw common::NotImplemented( FromHere(), "mononomial_exponents not implemented");
  static RealMatrix dummy;
  return dummy;
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3
