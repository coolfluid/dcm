// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_navierstokes_BCWall_hpp
#define cf3_dcm_equations_navierstokes_BCWall_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/navierstokes/LibNavierStokes.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

/// @brief Viscous adiabatic moving wall boundary condition for Navier-Stokes equations
///
/// The velocity inside the wall is reversed, so that a velocity flux results in zero.
/// The wall temperature is interpolated from inside.
///
/// By default the wall is not moving.
class dcm_equations_navierstokes_API BCWall1D : public dcm::core::BC<1,3>
{
public:
  BCWall1D(const std::string& name) : dcm::core::BC<1,3>(name) {}
  virtual ~BCWall1D() {}
  static std::string type_name() { return "BCWall1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution, 
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution );
};

////////////////////////////////////////////////////////////////////////////////

/// @brief Viscous adiabatic moving wall boundary condition for Navier-Stokes equations
///
/// The velocity inside the wall is reversed, so that a velocity flux results in zero.
/// The wall temperature is interpolated from inside.
///
/// By default the wall is not moving.
class dcm_equations_navierstokes_API BCWall2D : public dcm::core::BC<2,4>
{
public:
  BCWall2D(const std::string& name);
  virtual ~BCWall2D() {}
  static std::string type_name() { return "BCWall2D"; }
  
   virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                           const Matrix_NDIMxNEQS& inner_solution_gradient,
                                           const ColVector_NDIM& coords,
                                           const ColVector_NDIM& face_normal,
                                           RowVector_NEQS& boundary_solution,
                                           Matrix_NDIMxNEQS& boundary_solution_gradient );
private:
  Real m_wall_velocity;
  Real m_gamma;
  Real m_R;
};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_navierstokes_BCWall_hpp
