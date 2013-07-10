// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the CellConnectivitys of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "cf3/common/List.hpp"
#include "cf3/common/Table.hpp"
#include "cf3/common/Log.hpp"
#include "cf3/mesh/FaceCellConnectivity.hpp"

#include "cf3/dcm/core/CellConnectivity.hpp"
#include "cf3/mesh/ElementConnectivity.hpp"
#include "cf3/mesh/ElementType.hpp"
#include "cf3/mesh/Cells.hpp"

using namespace cf3::common;
using namespace cf3::mesh;

namespace cf3 {
namespace dcm {
namespace core {

/////////////////////////////////////////////////////////////////////////////////////

CellConnectivity::CellConnectivity ( const std::string& name ) :
  cf3::common::Component(name)
{
}

/////////////////////////////////////////////////////////////////////////////////////

CellConnectivity::~CellConnectivity()
{
}

/////////////////////////////////////////////////////////////////////////////////////

void CellConnectivity::compute_cell( const mesh::Cells& cells, const Uint elem_idx)
{
  if (m_cell.comp != &cells || m_cell.idx != elem_idx)
  {
    m_cell = Entity(cells,elem_idx);
    cf3_assert( is_not_null( m_cell.comp->connectivity_cell2face() ) );
    const ElementConnectivity& face_connectivity = *m_cell.comp->connectivity_cell2face();
    const Uint nb_faces = m_cell.element_type().nb_faces();
    m_faces.resize(nb_faces);
    m_orientations.resize(nb_faces);
    m_rotations.resize(nb_faces);
    m_is_bdry_face.resize(nb_faces);
    m_neighbor_cells.resize(nb_faces);
    m_neighbor_orientations.resize(nb_faces);
    m_neighbor_rotations.resize(nb_faces);
    m_neighbor_face_nb.resize(nb_faces);
    for (Uint face_nb=0; face_nb<nb_faces; ++face_nb)
    {
      cf3_assert( m_cell.idx < face_connectivity.size() );
      cf3_assert( face_nb < face_connectivity.row_size() );
      m_faces[face_nb] = face_connectivity[m_cell.idx][face_nb];
      Entity& face = m_faces[face_nb];
      cf3_assert( is_not_null( face.comp->connectivity_face2cell() ) );
      const FaceCellConnectivity& cell_connectivity = *face.comp->connectivity_face2cell();
      cf3_assert( face.idx < cell_connectivity.is_bdry_face().size() );
      if (cell_connectivity.is_bdry_face()[face.idx])
      {
        if ( m_cell != cell_connectivity.connectivity()[face.idx][LEFT] )
          throw InvalidStructure(FromHere(), "Face to Cell connectivity is wrong for boundary");

        m_is_bdry_face[face_nb] = true;

        m_neighbor_cells[face_nb] = face;
        m_neighbor_face_nb[face_nb] = 0;

        m_rotations[face_nb] = cell_connectivity.cell_rotation()[face.idx][LEFT];
        m_orientations[face_nb] = cell_connectivity.cell_orientation()[face.idx][LEFT];

        m_neighbor_orientations[face_nb] = MATCHED;
        m_neighbor_rotations[face_nb] = 0;
      }
      else
      {
        m_is_bdry_face[face_nb] = false;

        cf3_assert(face.idx < cell_connectivity.connectivity().size());
        cf3_assert(is_not_null(cell_connectivity.connectivity()[face.idx][LEFT].comp));
        cf3_assert(is_not_null(cell_connectivity.connectivity()[face.idx][RIGHT].comp));
        if ( m_cell == cell_connectivity.connectivity()[face.idx][LEFT] )
        {
          m_orientations[face_nb] = MATCHED;
          m_neighbor_orientations[face_nb] = INVERTED;
          m_rotations[face_nb] = cell_connectivity.cell_rotation()[face.idx][LEFT];
          m_neighbor_rotations[face_nb] = cell_connectivity.cell_rotation()[face.idx][RIGHT];
          m_neighbor_cells[face_nb] = cell_connectivity.connectivity()[face.idx][RIGHT];
          m_neighbor_face_nb[face_nb] = cell_connectivity.face_number()[face.idx][RIGHT];
        }
        else
        {
          m_orientations[face_nb] = INVERTED;
          m_neighbor_orientations[face_nb] = MATCHED;
          m_rotations[face_nb] = cell_connectivity.cell_rotation()[face.idx][RIGHT];
          m_neighbor_rotations[face_nb] = cell_connectivity.cell_rotation()[face.idx][LEFT];;
          m_neighbor_cells[face_nb] = cell_connectivity.connectivity()[face.idx][LEFT];
          m_neighbor_face_nb[face_nb] = cell_connectivity.face_number()[face.idx][LEFT];
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3

/////////////////////////////////////////////////////////////////////////////////////
