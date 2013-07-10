// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_core_CellConnectivity_hpp
#define cf3_dcm_core_CellConnectivity_hpp

#include "cf3/mesh/Entities.hpp"
#include "cf3/dcm/core/LibCore.hpp"

namespace cf3 {
namespace mesh { class Cells; }
namespace dcm {
namespace core {

/////////////////////////////////////////////////////////////////////////////////////

class dcm_core_API CellConnectivity : public cf3::common::Component {

public: // types

  enum FaceOrientation {MATCHED=0, INVERTED=1};

public: // functions

  /// Contructor
  /// @param name of the component
  CellConnectivity ( const std::string& name );

  /// Virtual destructor
  virtual ~CellConnectivity();

  /// Get the class name
  static std::string type_name () { return "CellConnectivity"; }

  void compute_cell( const mesh::Cells& cells, const Uint elem_idx);

  const mesh::Entity& cell() const { return m_cell; }

  const std::vector<Uint>& rotations() const { return m_rotations; }

  const std::vector<Uint>& orientations() const { return m_orientations; }
  
  const std::vector<mesh::Entity>& neighbor_cells() const { return m_neighbor_cells; }

  const std::vector<Uint>& neighbour_rotations() const { return m_neighbor_rotations; }

  const std::vector<Uint>& neighbour_orientations() const { return m_neighbor_orientations; }

  const std::vector<Uint>& neighbour_face_nb() const { return m_neighbor_face_nb; }

  const std::vector<bool>& is_bdry_face() const { return m_is_bdry_face; }

private:
  
  mesh::Entity              m_cell;
  std::vector<Uint>         m_orientations;
  std::vector<Uint>         m_rotations;
  std::vector<mesh::Entity> m_faces;
  std::vector<bool>         m_is_bdry_face;
  std::vector<mesh::Entity> m_neighbor_cells;
  std::vector<Uint>         m_neighbor_orientations;
  std::vector<Uint>         m_neighbor_rotations;
  std::vector<Uint>         m_neighbor_face_nb;
};

/////////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3

#endif // cf3_dcm_core_CellConnectivity_hpp
