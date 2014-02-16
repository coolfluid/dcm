// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file BCSubsonicInlet2D.hpp
/// @brief Two Navierstokes/Euler subsonic inlet boundary conditions
///
/// - BCSubsonicInlet2D
/// - BCSubsonicInletUT2D

#ifndef cf3_dcm_equations_euler_BCSubsonicInlet2D_hpp
#define cf3_dcm_equations_euler_BCSubsonicInlet2D_hpp

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
/// Default configuration: Tt = 25 Celsius, Pt = 1 bar , angle = 0 rad, gas = air
class dcm_equations_euler_API BCSubsonicInlet2D : public dcm::core::BC<2,4>
{
private:
  enum {Rho=0, RhoUx=1, RhoUy=2, RhoE=3};
  
public:
  BCSubsonicInlet2D(const std::string& name);
  virtual ~BCSubsonicInlet2D() {}
  static std::string type_name() { return "BCSubsonicInlet2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
  
private:
  void config_Tt();
  void config_Pt();
  void config_alpha();
  
private:
  math::AnalyticalFunction m_function_Tt;
  math::AnalyticalFunction m_function_Pt;
  math::AnalyticalFunction m_function_alpha;

  Real m_Tt;
  Real m_Pt;
  Real m_alpha;

  Real m_R;
  Real m_gamma;
  ColVector_NDIM m_U;

  Real m_T_inner;
  Real m_rho_inner;
  Real m_p_inner;
  Real m_rhoE_inner;
  Real m_uuvv_inner;
  Real m_M2_inner;
  Real m_coeff_inner;
  Real m_pow_coeff_inner;

  Real m_M;
  Real m_tan_alpha;
  Real m_T;
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

#endif // cf3_dcm_BCSubsonicInlet2D_hpp
