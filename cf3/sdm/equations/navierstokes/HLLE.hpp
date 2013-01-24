// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_navierstokes_HLLE_hpp
#define cf3_sdm_equations_navierstokes_HLLE_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/sdm/equations/navierstokes/Terms1D.hpp"
#include "cf3/sdm/equations/navierstokes/Terms2D.hpp"
#include "cf3/solver/RiemannSolver.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_navierstokes_API HLLE1D : public solver::RiemannSolver<Terms1D>
{
public:

  HLLE1D(const std::string& name) : solver::RiemannSolver<Terms1D>(name) {}
  virtual ~HLLE1D() {}
  static std::string type_name() { return "HLLE1D"; }

  virtual void compute_riemann_flux( const Terms1D::DATA& left, const Terms1D::DATA& right, const Terms1D::ColVector_NDIM& normal,
                                     Terms1D::RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler1d::compute_hlle_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_navierstokes_API HLLE2D : public solver::RiemannSolver<Terms2D>
{
public:

  HLLE2D(const std::string& name) : solver::RiemannSolver<Terms2D>(name) {}
  virtual ~HLLE2D() {}
  static std::string type_name() { return "HLLE2D"; }

  virtual void compute_riemann_flux( const Terms2D::DATA& left, const Terms2D::DATA& right, const Terms2D::ColVector_NDIM& normal,
                                     Terms2D::RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_hlle_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_navierstokes_HLLE_hpp
