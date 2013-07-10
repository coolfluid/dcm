// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file BCPressureVelocityInlet2D.hpp
/// @brief Two Navierstokes/Euler subsonic inlet boundary conditions
///
/// - BCPressureVelocityInlet2D
/// - BCPressureVelocityInletUT2D

#ifndef cf3_dcm_equations_euler_BCPressureVelocityInlet2D_hpp
#define cf3_dcm_equations_euler_BCPressureVelocityInlet2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/math/AnalyticalFunction.hpp"
#include "cf3/math/Functions.hpp"

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/euler/LibEuler.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

/// @brief Subsonic inlet given total temperature, total pressure and flow angle
/// @author Willem Deconinck
/// Default configuration: T = 25 Celsius, Pt = 1 bar , angle = 0 rad, gas = air
class dcm_equations_euler_API BCPressureVelocityInlet2D : public dcm::core::BC<2,4>
{
private:
  enum {Rho=0, RhoUx=1, RhoUy=2, RhoE=3};
  
public:
  BCPressureVelocityInlet2D(const std::string& name);
  virtual ~BCPressureVelocityInlet2D() {}
  static std::string type_name() { return "BCPressureVelocityInlet2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
  
private:
  void config_p();
  void config_u();
  
private:
  math::AnalyticalFunction m_function_p;
  math::AnalyticalFunction m_function_u;

  Real m_R;
  Real m_gamma;
  ColVector_NDIM m_U;
  Real m_p;
  Real m_rho;
  Real m_uuvv;
  Real m_rhoE;
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_BCPressureVelocityInlet2D_hpp
