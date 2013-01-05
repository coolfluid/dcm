// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_core_FaceConnectivity_hpp
#define cf3_sdm_core_FaceConnectivity_hpp

#include "cf3/mesh/Entities.hpp"
#include "cf3/sdm/core/LibCore.hpp"

namespace cf3 {
namespace mesh { class Cells; class Faces; }
namespace sdm {
namespace core {

/////////////////////////////////////////////////////////////////////////////////////

class sdm_core_API FaceConnectivity : public cf3::common::Component {

public: // types

  enum FaceOrientation {MATCHED=0, INVERTED=1};

public: // functions

  /// Contructor
  /// @param name of the component
  FaceConnectivity ( const std::string& name );

  /// Virtual destructor
  virtual ~FaceConnectivity();

  /// Get the class name
  static std::string type_name () { return "FaceConnectivity"; }

  void compute_face( const mesh::Faces& faces, const Uint face_idx);

  const mesh::Entity& face() const { return m_face; }
  
  bool is_bdry_face() const { return m_is_bdry_face; }
  
  Uint rotation() const { return m_rotation; }

  Uint orientation() const { return m_orientation; }
  
  const std::vector<mesh::Entity>& cells() const { return m_cells; }
  
  const std::vector<Uint>& cells_rotation() const { return m_cells_rotation; }

  const std::vector<Uint>& cells_orientation() const { return m_cells_orientation; }

  const std::vector<Uint>& cells_face_nb() const { return m_cells_face_nb; }

private:
  
  mesh::Entity                m_face;
  bool                        m_is_bdry_face;
  Uint                        m_orientation;
  Uint                        m_rotation;
  std::vector<mesh::Entity>   m_cells;
  std::vector<Uint>           m_cells_orientation;
  std::vector<Uint>           m_cells_rotation;
  std::vector<Uint>           m_cells_face_nb;
};

/////////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

#endif // cf3_sdm_core_FaceConnectivity_hpp
