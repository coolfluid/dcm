// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/sdm/equations/navierstokes/BCWall.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCWall1D,solver::BC,LibNavierStokes> BCWall1D_builder;
common::ComponentBuilder<BCWall2D,solver::BC,LibNavierStokes> BCWall2D_builder;

////////////////////////////////////////////////////////////////////////////////

void BCWall1D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
{
  boundary_solution[0] =  inner_solution[0];
  boundary_solution[1] = -inner_solution[1];
  boundary_solution[2] =  inner_solution[2];
}

////////////////////////////////////////////////////////////////////////////////

BCWall2D::BCWall2D(const std::string& name) :
  sdm::core::BC<2,4>(name),
  m_wall_velocity(0.)
{
  options().add("wall_velocity",m_wall_velocity)
      .description("The velocity of the wall")
      .link_to(&m_wall_velocity)
      .mark_basic();
}
void BCWall2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
{
  boundary_solution[0]=  inner_solution[0];
  boundary_solution[1]= -inner_solution[1] - 2.*inner_solution[0]*face_normal[YY]*m_wall_velocity;
  boundary_solution[2]= -inner_solution[2] + 2.*inner_solution[0]*face_normal[XX]*m_wall_velocity;
  boundary_solution[3]=  inner_solution[3];
}



/* FVM 3D  in: Vector3D& wall_velocity
 * ------
    double  Wall_velocity_tang ;
    Vector3D ur_norm, ur_tang, vr_tot, uw_tang;

    ur_norm = dot(Win.v, norm_dir)*norm_dir;
    ur_tang = Win.v - ur_norm;

    uw_tang = wall_velocity - dot(norm_dir,wall_velocity)*norm_dir;

    ur_norm = -ur_norm;
    ur_tang = 2.0*uw_tang - ur_tang;
    vr_tot = ur_norm + ur_tang;

    Temp.v = vr_tot;

    //  Fixed Wall Temperature or constant extrapolation for Adiabatic
    if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
        if (pressure_gradient != Vector3D_ZERO){
            Temp.rho = Wout.p/(Temp.R*Wout.T());
        } else {
            Temp.rho = Win.p/(Temp.R*Wout.T());
    }
*/

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // sdm
} // cf3
