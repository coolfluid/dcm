// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_euler_Terms_hpp
#define cf3_dcm_equations_euler_Terms_hpp

////////////////////////////////////////////////////////////////////////////////

#include <boost/mpl/vector.hpp>

#include "cf3/dcm/equations/euler/RightHandSide1D.hpp"
#include "cf3/dcm/equations/euler/RightHandSide2D.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

typedef boost::mpl::vector<
  RightHandSide1D,
  RightHandSide2D
> Terms;

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_euler_Terms_hpp
