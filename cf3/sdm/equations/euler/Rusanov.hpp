// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_euler_Rusanov_hpp
#define cf3_sdm_equations_euler_Rusanov_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/sdm/equations/euler/Terms1D.hpp"
#include "cf3/sdm/equations/euler/Terms2D.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API Rusanov1D : public solver::RiemannSolver<Terms1D>
{
public:

  Rusanov1D(const std::string& name) : solver::RiemannSolver<Terms1D>(name) {}
  virtual ~Rusanov1D() {}
  static std::string type_name() { return "Rusanov1D"; }

  virtual void compute_riemann_flux( const Terms1D::DATA& left, const Terms1D::DATA& right, const Terms1D::ColVector_NDIM& normal,
                                     Terms1D::RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler1d::compute_rusanov_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API Rusanov2D : public solver::RiemannSolver<Terms2D>
{
public:

  Rusanov2D(const std::string& name) : solver::RiemannSolver<Terms2D>(name) {}
  virtual ~Rusanov2D() {}
  static std::string type_name() { return "Rusanov2D"; }

  virtual void compute_riemann_flux( const Terms2D::DATA& left, const Terms2D::DATA& right, const Terms2D::ColVector_NDIM& normal,
                                     Terms2D::RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_rusanov_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_euler_Rusanov_hpp
