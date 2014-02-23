// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_les_ComputeSubFilterScale_hpp
#define cf3_dcm_equations_les_ComputeSubFilterScale_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/common/Action.hpp"
#include "cf3/dcm/equations/les/LibLES.hpp"

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
namespace cf3 {
  namespace mesh {
    class Field;
    class Region;
  }
  namespace dcm {
    namespace equations {
      namespace les {
        class EddyViscosityModel;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

/// @brief Compute Sub-Filter-Scale values
///
/// This component computes the modeled values for
///   - sub-filter-scale kinetic energy
///   - sub-filter-scale kinematic viscosity a.k.a. eddy viscosity
///   - sub-filter-scale heat conduction
///
/// This component needs
///   - density
///   - velocity
/// @author Willem Deconinck

class dcm_equations_les_API ComputeSubFilterScale : public common::Action {

public: // functions

  /// Contructor
  /// @param name of the component
  ComputeSubFilterScale ( const std::string& name );

  /// Virtual destructor
  virtual ~ComputeSubFilterScale();

  /// Get the class name
  static std::string type_name () { return "ComputeSubFilterScale"; }

  virtual void execute();

private:

  Handle<mesh::Field> m_density;
  Handle<mesh::Field> m_velocity;
  Handle<mesh::Field> m_k_sfs;
  Handle<mesh::Field> m_kappaT;
  Handle<mesh::Field> m_nuT;
  
  Handle< EddyViscosityModel > m_sfs_model;
};

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_les_ComputeSubFilterScale_hpp
