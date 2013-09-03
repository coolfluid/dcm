// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_euler_RiemannSolvers_hpp
#define cf3_dcm_equations_euler_RiemannSolvers_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/physics/euler/euler1d/Data.hpp"
#include "cf3/physics/euler/euler1d/Functions.hpp"
#include "cf3/physics/euler/euler2d/Data.hpp"
#include "cf3/physics/euler/euler2d/Functions.hpp"
#include "cf3/dcm/equations/euler/LibEuler.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

typedef solver::RiemannSolver<
    cf3::physics::euler::euler1d::Data,
    cf3::physics::euler::euler1d::NDIM,
    cf3::physics::euler::euler1d::NEQS > RiemannSolver1D;

typedef solver::RiemannSolver<
    cf3::physics::euler::euler2d::Data,
    cf3::physics::euler::euler2d::NDIM,
    cf3::physics::euler::euler2d::NEQS > RiemannSolver2D;

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API Roe1D : public RiemannSolver1D
{
public:

  Roe1D(const std::string& name) : RiemannSolver1D(name) {}
  virtual ~Roe1D() {}
  static std::string type_name() { return "Roe1D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler1d::compute_roe_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API Roe2D : public RiemannSolver2D
{
public:

  Roe2D(const std::string& name) : RiemannSolver2D(name) {}
  virtual ~Roe2D() {}
  static std::string type_name() { return "Roe2D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_roe_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API HLLE1D : public RiemannSolver1D
{
public:

  HLLE1D(const std::string& name) : RiemannSolver1D(name) {}
  virtual ~HLLE1D() {}
  static std::string type_name() { return "HLLE1D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler1d::compute_hlle_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API HLLE2D : public RiemannSolver2D
{
public:

  HLLE2D(const std::string& name) : RiemannSolver2D(name) {}
  virtual ~HLLE2D() {}
  static std::string type_name() { return "HLLE2D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_hlle_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API Rusanov1D : public RiemannSolver1D
{
public:

  Rusanov1D(const std::string& name) : RiemannSolver1D(name) {}
  virtual ~Rusanov1D() {}
  static std::string type_name() { return "Rusanov1D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler1d::compute_rusanov_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API Rusanov2D : public RiemannSolver2D
{
public:

  Rusanov2D(const std::string& name) : RiemannSolver2D(name) {}
  virtual ~Rusanov2D() {}
  static std::string type_name() { return "Rusanov2D"; }

  virtual void compute_riemann_flux( const Data& left, const Data& right, const ColVector_NDIM& normal,
                                     RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_rusanov_flux(left,right,normal,flux,wave_speed);
  }
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_euler_RiemannSolvers_hpp
