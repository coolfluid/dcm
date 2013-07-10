// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/dcm/equations/euler/BCMirror.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<BCMirror1D,solver::BC,LibEuler> BCMirror1D_builder;
common::ComponentBuilder<BCMirror2D,solver::BC,LibEuler> BCMirror2D_builder;
common::ComponentBuilder<BCMirror3D,solver::BC,LibEuler> BCMirror3D_builder;

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
  static ColVector_NDIM rhoU;
    
  // velocity on the inside of the face
  rhoU[XX] = inner_solution[1];
  rhoU[YY] = inner_solution[2];

  // velocity in outward_normal direction
  static Real rhoU_normal;
  rhoU_normal = rhoU.transpose()*face_normal;

  // Modify velocity to become the outside velocity of the face,
  // and being the mirror of the inside velocity
  rhoU.noalias() -= 2.*rhoU_normal*face_normal;

  // Set the outside boundary state
  boundary_solution[0]=inner_solution[0];
  boundary_solution[1]=rhoU[XX];
  boundary_solution[2]=rhoU[YY];
  boundary_solution[3]=inner_solution[3];
}

////////////////////////////////////////////////////////////////////////////////

void BCMirror3D::compute_boundary_solution( const RowVector_NEQS& inner_solution, 
                                            const ColVector_NDIM& coords,
                                            const ColVector_NDIM& face_normal, 
                                            RowVector_NEQS& boundary_solution )
{
  static ColVector_NDIM rhoU;
    
  // velocity on the inside of the face
  rhoU[XX] = inner_solution[1];
  rhoU[YY] = inner_solution[2];
  rhoU[ZZ] = inner_solution[3];

  // velocity in outward_normal direction
  static Real rhoU_normal;
  rhoU_normal = rhoU.transpose()*face_normal;

  // Modify velocity to become the outside velocity of the face,
  // and being the mirror of the inside velocity
  rhoU.noalias() -= 2.*rhoU_normal*face_normal;

  // Set the outside boundary state
  boundary_solution[0]=inner_solution[0];
  boundary_solution[1]=rhoU[XX];
  boundary_solution[2]=rhoU[YY];
  boundary_solution[3]=rhoU[ZZ];
  boundary_solution[4]=inner_solution[4];
}

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3
