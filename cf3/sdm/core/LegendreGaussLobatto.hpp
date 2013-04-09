// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file LegendreGaussLobatto.hpp
///
/// This file defines the distribution for solution points and flux points
/// for a 1D shape function, used in tensorial elements.
/// The flux points are located at the Legendre-Gauss-Lobatto flux-points,
/// and the solution points coincide as much as possible.

#ifndef cf3_sdm_core_LegendreGaussLobatto_hpp
#define cf3_sdm_core_LegendreGaussLobatto_hpp


#include "cf3/math/MatrixTypes.hpp"
#include "cf3/sdm/core/LibCore.hpp"

#include <boost/assign/list_of.hpp>
using namespace boost::assign;

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

struct sdm_core_API LegendreGaussLobatto
{
  LegendreGaussLobatto(const Uint p)
  {
    nb_sol_pts = p+1;
    nb_flx_pts = p+2;
    
    sol_pts.resize(nb_sol_pts);
    flx_pts.resize(nb_flx_pts);

    switch (p)
    {
    case 0:
      flx_pts << -1,1;
      break;
    case 1:
      flx_pts << -1,0,1;
      break;
    case 2:
      flx_pts << -1,-1/sqrt(3.),1/sqrt(3.),1;
      break;
    case 3:
      flx_pts << -1,-sqrt(3./5.),0,+sqrt(3./5.),1;
      break;
    case 4:
      flx_pts << -1, -sqrt((3.+2.*sqrt(6./5.))/7.), -sqrt((3.-2.*sqrt(6./5.))/7.), +sqrt((3.-2.*sqrt(6./5.))/7.), +sqrt((3.+2.*sqrt(6./5.))/7.), 1;
      break;
    case 5:
      flx_pts << -1, -sqrt(5.+2.*sqrt(10./7.))/3., -sqrt(5.-2.*sqrt(10./7.))/3., 0., +sqrt(5.-2.*sqrt(10./7.))/3., +sqrt(5.+2.*sqrt(10./7.))/3., +1;
      break;
    default:
      throw common::NotImplemented(FromHere(),"1D flux-point locations for P"+common::to_str(p)+" are not yet defined");
      break;
    }

    // Make solution points symmetric and collocated as much as possible with flux points
    Uint s = 0;
    for (Uint f = 0; f < sol_pts.size()/2; ++f, ++s)
      sol_pts[s] = flx_pts[f];
    if (sol_pts.size()%2 == 1)
    {
      sol_pts[s] = 0.5*(flx_pts[sol_pts.size()/2]+flx_pts[sol_pts.size()/2+1]);
      ++s;
    }
    for (Uint f = flx_pts.size()/2+1; f < flx_pts.size(); ++f, ++s)
      sol_pts[s] = flx_pts[f];
  }

  RealVector sol_pts;
  RealVector flx_pts;
  Uint nb_sol_pts;
  Uint nb_flx_pts;
};

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

#endif // cf3_sdm_core_LegendreGaussLobatto_hpp
