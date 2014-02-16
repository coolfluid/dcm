// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_Terms_hpp
#define cf3_dcm_equations_lineuler_Terms_hpp

////////////////////////////////////////////////////////////////////////////////

#include <boost/mpl/vector.hpp>

#include "cf3/dcm/equations/lineuler/RightHandSide2D.hpp"
#include "cf3/dcm/equations/lineuler/SourceMonopoleUniform.hpp"
#include "cf3/dcm/equations/lineuler/SourceDipole.hpp"
#include "cf3/dcm/equations/lineuler/SourceQuadrupole.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

typedef boost::mpl::vector<
  RightHandSide2D,
  SourceMonopoleUniform2D,
  SourceMonopoleUniform3D,
  SourceDipole2D,
  SourceQuadrupole2D
> Terms;

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_Terms_hpp
