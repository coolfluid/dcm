// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_navierstokes_ComputeBoundaryLayer_hpp
#define cf3_dcm_equations_navierstokes_ComputeBoundaryLayer_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/equations/navierstokes/LibNavierStokes.hpp"
#include "cf3/common/Action.hpp"

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
namespace cf3 {
  namespace mesh {
    class Field;
    class Region;
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

/// @brief Compute Boundary layer values
///
/// This component computes the values for
///   - y-plus
///   - tau_w (wall shear stress)
///   - ustar (friction velocity)
///
/// This component needs
///   - y0: a value of the wall-distance of the first cell
///   - density: The density field, defined in the domain
///   - velocity: The velocity field, defined in the domain
///
/// The values are only computed for points in the faces
/// of the wall.
/// Fields for yplus, tau_w and ustar need to be created
/// first in the faces of the wall, and configured in this
/// component
///
/// @todo 3D implementation
///
/// @author Willem Deconinck

class dcm_equations_navierstokes_API ComputeBoundaryLayer : public common::Action {

public: // functions

  /// Contructor
  /// @param name of the component
  ComputeBoundaryLayer ( const std::string& name );

  /// Virtual destructor
  virtual ~ComputeBoundaryLayer();

  /// Get the class name
  static std::string type_name () { return "ComputeBoundaryLayer"; }

  virtual void execute();

private:
  
  std::vector< Handle<mesh::Region> > m_wall_regions;
  Handle<mesh::Field> m_velocity;
  Handle<mesh::Field> m_density;
  Handle<mesh::Field> m_yplus;
  Handle<mesh::Field> m_tau;
  Handle<mesh::Field> m_ustar;
};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_navierstokes_ComputeBoundaryLayer_hpp
