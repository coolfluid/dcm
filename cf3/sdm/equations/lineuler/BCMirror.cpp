// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/sdm/equations/lineuler/BCMirror.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCMirror1D,solver::BC,LibLinEuler> BCMirror1D_builder;
common::ComponentBuilder<BCMirror2D,solver::BC,LibLinEuler> BCMirror2D_builder;
common::ComponentBuilder<BCMirror3D,solver::BC,LibLinEuler> BCMirror3D_builder;

////////////////////////////////////////////////////////////////////////////////

void BCMirror1D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                            const ColVector_NDIM& coords,
                                            const ColVector_NDIM& face_normal,
                                            RowVector_NEQS& boundary_solution )
{
  boundary_solution[0] =   inner_solution[0];
  boundary_solution[1] = - inner_solution[1];
  boundary_solution[2] =   inner_solution[2];
}

////////////////////////////////////////////////////////////////////////////////

void BCMirror2D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                            const ColVector_NDIM& coords,
                                            const ColVector_NDIM& face_normal, 
                                            RowVector_NEQS& boundary_solution )
{
  static ColVector_NDIM rho0_U;
    
  // velocity on the inside of the face
  rho0_U[XX] = inner_solution[1];
  rho0_U[YY] = inner_solution[2];

  // velocity in outward_normal direction
  static Real rho0_U_normal;
  rho0_U_normal = rho0_U.transpose()*face_normal;

  // Modify velocity to become the outside velocity of the face,
  // and being the mirror of the inside velocity
  rho0_U.noalias() -= 2.*rho0_U_normal*face_normal;

  // Set the outside boundary state
  boundary_solution[0]=inner_solution[0];
  boundary_solution[1]=rho0_U[XX];
  boundary_solution[2]=rho0_U[YY];
  boundary_solution[3]=inner_solution[3];
}

////////////////////////////////////////////////////////////////////////////////

void BCMirror3D::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                            const ColVector_NDIM& coords,
                                            const ColVector_NDIM& face_normal, 
                                            RowVector_NEQS& boundary_solution )
{
  static ColVector_NDIM rho0_U;
    
  // velocity on the inside of the face
  rho0_U[XX] = inner_solution[1];
  rho0_U[YY] = inner_solution[2];
  rho0_U[ZZ] = inner_solution[3];

  // velocity in outward_normal direction
  static Real rho0_U_normal;
  rho0_U_normal = rho0_U.transpose()*face_normal;

  // Modify velocity to become the outside velocity of the face,
  // and being the mirror of the inside velocity
  rho0_U.noalias() -= 2.*rho0_U_normal*face_normal;

  // Set the outside boundary state
  boundary_solution[0]=inner_solution[0];
  boundary_solution[1]=rho0_U[XX];
  boundary_solution[2]=rho0_U[YY];
  boundary_solution[3]=rho0_U[ZZ];
  boundary_solution[4]=inner_solution[4];
}

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
