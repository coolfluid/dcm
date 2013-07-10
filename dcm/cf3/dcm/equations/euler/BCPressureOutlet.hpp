// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file BCPressureOutlet2D.hpp
/// @brief Two Navierstokes/Euler subsonic inlet boundary conditions
///
/// - BCPressureOutlet2D
/// - BCPressureOutletUT2D

#ifndef cf3_dcm_equations_euler_BCPressureOutlet2D_hpp
#define cf3_dcm_equations_euler_BCPressureOutlet2D_hpp

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

/// @brief Subsonic outlet given a pressure
/// @author Willem Deconinck
class dcm_equations_euler_API BCPressureOutlet2D : public dcm::core::BC<2,4>
{
private:
  enum {Rho=0, RhoUx=1, RhoUy=2, RhoE=3};
  
public:
  BCPressureOutlet2D(const std::string& name);
  virtual ~BCPressureOutlet2D() {}
  static std::string type_name() { return "BCPressureOutlet2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
   
private:
  void config_p();
  
private:
  math::AnalyticalFunction m_function_p;

  Real m_p;
  Real m_gamma;

  Real m_u_inner;
  Real m_v_inner;
  Real m_rho_inner;
  Real m_uuvv_inner;

};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_BCPressureOutlet2D_hpp
