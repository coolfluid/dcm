// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the FaceConnectivitys of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "cf3/common/List.hpp"
#include "cf3/common/Table.hpp"
#include "cf3/mesh/FaceCellConnectivity.hpp"
#include "cf3/mesh/ElementConnectivity.hpp"
#include "cf3/mesh/ElementType.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Faces.hpp"
#include "cf3/sdm/core/FaceConnectivity.hpp"

using namespace cf3::common;
using namespace cf3::mesh;

namespace cf3 {
namespace sdm {
namespace core {

/////////////////////////////////////////////////////////////////////////////////////

FaceConnectivity::FaceConnectivity ( const std::string& name ) :
  cf3::common::Component(name)
{
  m_cells.resize(2);
  m_cells_rotation.resize(2);
  m_cells_orientation.resize(2);
  m_cells_face_nb.resize(2);
}

/////////////////////////////////////////////////////////////////////////////////////

FaceConnectivity::~FaceConnectivity()
{
}

/////////////////////////////////////////////////////////////////////////////////////

void FaceConnectivity::compute_face( const mesh::Faces& faces, const Uint face_idx)
{
  if (m_face.comp != &faces || m_face.idx != face_idx)
  {
    m_face = Entity(faces,face_idx);
    cf3_assert( is_not_null( m_face.comp->connectivity_face2cell() ) );
    const mesh::FaceCellConnectivity& cell_connectivity = *m_face.comp->connectivity_face2cell();
    cf3_assert( m_face.idx < cell_connectivity.size() );

    m_orientation             = MATCHED;
    m_rotation                = 0;
    m_cells[LEFT]             = cell_connectivity.connectivity()[m_face.idx][LEFT];
    m_cells_rotation[LEFT]    = cell_connectivity.cell_rotation()[m_face.idx][LEFT];
    m_cells_orientation[LEFT] = cell_connectivity.cell_orientation()[m_face.idx][LEFT];
    m_cells_face_nb[LEFT]     = cell_connectivity.face_number()[m_face.idx][LEFT];
    
    m_is_bdry_face = cell_connectivity.is_bdry_face()[m_face.idx];
    if (m_is_bdry_face == false)
    {
      m_cells[RIGHT]             = cell_connectivity.connectivity()[m_face.idx][RIGHT];
      m_cells_rotation[RIGHT]    = cell_connectivity.cell_rotation()[m_face.idx][RIGHT];
      m_cells_orientation[RIGHT] = cell_connectivity.cell_orientation()[m_face.idx][RIGHT];
      m_cells_face_nb[RIGHT]     = cell_connectivity.face_number()[m_face.idx][RIGHT];
    } 
  }
}

/////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

/////////////////////////////////////////////////////////////////////////////////////
