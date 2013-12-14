// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_core_BC_hpp
#define cf3_dcm_core_BC_hpp

#include "cf3/physics/MatrixTypes.hpp"

#include <cf3/common/Log.hpp>

#include "cf3/math/Consts.hpp"

#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Faces.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Connectivity.hpp"

#include "cf3/solver/BC.hpp"

#include "cf3/dcm/core/LibCore.hpp"
#include "cf3/dcm/core/FaceConnectivity.hpp"
#include "cf3/dcm/core/ShapeFunction.hpp"
#include "cf3/dcm/core/Reconstructions.hpp"
#include "cf3/dcm/core/Metrics.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM, Uint NB_EQS>
class dcm_core_API BC : public cf3::solver::BC {
public: // types

  enum {NDIM  = NB_DIM};
  enum {NEQS  = NB_EQS};

  typedef typename physics::MatrixTypes<NDIM,NEQS>::ColVector_NDIM ColVector_NDIM;
  typedef typename physics::MatrixTypes<NDIM,NEQS>::RowVector_NEQS RowVector_NEQS;
  typedef typename physics::MatrixTypes<NDIM,NEQS>::Matrix_NDIMxNEQS Matrix_NDIMxNEQS;
  typedef typename physics::MatrixTypes<NDIM,NDIM>::Matrix_NDIMxNDIM Matrix_NDIMxNDIM;

public: // functions
  /// Contructor
  /// @param name of the component
  BC ( const std::string& name );

  /// Virtual destructor
  virtual ~BC() {}

  /// Get the class name
  static std::string type_name () { return "BC"; }

  virtual void execute();
  
  virtual bool loop_faces(const Handle<mesh::Entities const>& cells);

  virtual void compute_element(const Uint elem_idx);

  virtual void compute_boundary_solution(const Uint cell_face_flx_pt, RowVector_NEQS& boundary_solution,
                                         Matrix_NDIMxNEQS& boundary_solution_gradient);

  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution, 
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution );

  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient );

  virtual void set_face_solution(const Uint face_pt, RowVector_NEQS& boundary_solution,
                                 Matrix_NDIMxNEQS& boundary_solution_gradient);
  
protected:

  Handle< dcm::core::FaceConnectivity >     m_connected;
  Handle< mesh::Faces const >               m_faces;
  Handle< mesh::Cells const >               m_cells;
  Handle< mesh::Space const >               m_faces_space;
  Handle< mesh::Space const >               m_cells_space;
  Handle< dcm::core::ShapeFunction const >  m_faces_sf;
  Handle< dcm::core::ShapeFunction const >  m_cells_sf;
  Handle< dcm::core::Metrics<NDIM> >        m_metrics;
  Handle< dcm::core::ElementMetrics<NDIM> > m_cell_metrics;

  Uint m_face_idx;
  Uint m_cell_idx;
  Uint m_cell_face_nb;
  Uint m_nb_face_pts;
  
  std::vector<RowVector_NEQS,Eigen::aligned_allocator<RowVector_NEQS> > m_boundary_solution;
  std::vector<Matrix_NDIMxNEQS,Eigen::aligned_allocator<Matrix_NDIMxNEQS> > m_boundary_solution_gradient;

  std::vector< std::vector<Uint> >  m_face_pts;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/////////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
BC<NB_DIM,NB_EQS>::BC ( const std::string& name ) :
  cf3::solver::BC(name)
{
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::execute()
{
  boost_foreach ( const Handle<Component>& region, m_regions )
  {
    if ( is_not_null(region) )
    {      
      boost_foreach( const mesh::Entities& faces, common::find_components_recursively< mesh::Entities >(*region) )
      {
        if ( loop_faces( faces.handle<mesh::Entities>() ) )
        {
          for (Uint face_idx=0; face_idx<faces.size(); ++face_idx)
            compute_element( face_idx );
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
bool BC<NB_DIM,NB_EQS>::loop_faces(const Handle<mesh::Entities const>& faces)
{
  if ( is_null(faces->handle<mesh::Faces>()) )
    return false;

  cf3_assert(m_fields);
  cf3_assert(m_bdry_fields);
  m_faces = faces->handle<mesh::Faces>();
  cf3_assert(m_faces);
  m_faces_space = m_bdry_fields->space(m_faces);
  cf3_assert(m_faces_space);
  m_faces_sf = m_faces_space->shape_function().handle<dcm::core::ShapeFunction>();
  cf3_assert(m_faces_sf);
  m_nb_face_pts = m_faces_sf->nb_sol_pts();
  m_boundary_solution.resize(m_nb_face_pts);
  m_boundary_solution_gradient.resize(m_nb_face_pts);
  // Create or find cellconnectivity
  if  ( Handle<Component const> found = m_faces_space->get_child("face_connectivity") )
    m_connected = const_cast<Component*>(found.get())
      ->handle< FaceConnectivity >();
  else
    m_connected = const_cast<mesh::Space*>(m_faces_space.get())
      ->create_component< FaceConnectivity >("face_connectivity");

  return true;
}
////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::compute_element( const Uint face_idx )
{
  m_face_idx = face_idx;
  // ---------------------------------------------------------------------------
  //                     COMPUTE CELL METRICS AND CONNECTIVITY
  // ---------------------------------------------------------------------------
  m_connected->compute_face(*m_faces,face_idx);
  cf3_assert(m_connected->is_bdry_face());
  m_cells = m_connected->cells()[LEFT].comp->handle<mesh::Cells>();
  cf3_assert(m_cells);
  m_cells_space = m_fields->space(m_cells);
  m_cell_idx = m_connected->cells()[LEFT].idx;
  cf3_assert(m_cells_space);
  m_cells_sf = m_cells_space->shape_function().handle<dcm::core::ShapeFunction>();
  cf3_assert(m_cells_sf);
  
  if  ( Handle<Component const> found = m_cells_space->get_child("metrics") )
  {
    m_metrics = const_cast<Component*>(found.get())->handle< dcm::core::Metrics<NDIM> >();
  }
  else
  {
    m_metrics = const_cast<mesh::Space*>(m_cells_space.get())->create_component< dcm::core::Metrics<NDIM> >("metrics");
    m_metrics->setup_for_space(m_cells_space);
  }
  m_cell_metrics = m_metrics->element(m_cell_idx);

  const std::vector<Uint>& cell_face_flx_pts = 
    m_cells_sf->face_flx_pts( m_connected->cells_face_nb()[LEFT],
                              m_connected->cells_orientation()[LEFT],
                              m_connected->cells_rotation()[LEFT] );
  
  // ---------------------------------------------------------------------------
  //                       COMPUTE BOUNDARY CONDITION
  // ---------------------------------------------------------------------------
  for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
  {
    const Uint cell_face_flx_pt = cell_face_flx_pts[face_pt];
      
    compute_boundary_solution(cell_face_flx_pt, m_boundary_solution[face_pt], m_boundary_solution_gradient[face_pt]);
    
    set_face_solution(face_pt,m_boundary_solution[face_pt], m_boundary_solution_gradient[face_pt]);
  }
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::compute_boundary_solution( const Uint cell_face_flx_pt, 
                                                   RowVector_NEQS& boundary_solution,
                                                   Matrix_NDIMxNEQS& boundary_solution_gradient )
{
  mesh::Connectivity::ConstRow nodes = m_cells_space->connectivity()[m_cell_idx];
  
  static RowVector_NEQS inner_solution;
  inner_solution.setZero();
  boost_foreach( const Uint sol_pt, m_metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).used_points() )
  {
    const Real C = m_metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).coeff(sol_pt);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      inner_solution[eq] += C * solution()->array()[nodes[sol_pt]][eq];
    }
  }
  static Matrix_NDIMxNEQS inner_solution_gradient;
  inner_solution_gradient.setZero();
  for( Uint d=0; d<NDIM; ++d )
  {
    boost_foreach( const Uint sol_pt, m_metrics->gradient_from_sol_pts_to_flx_pt(cell_face_flx_pt)[d].used_points() )
    {
      const Real C = m_metrics->gradient_from_sol_pts_to_flx_pt(cell_face_flx_pt)[d].coeff(sol_pt);
      for (Uint eq=0; eq<NEQS; ++eq)
      {
        inner_solution_gradient(d,eq) += C * solution()->array()[nodes[sol_pt]][eq];
      }
    }
  }
  inner_solution_gradient = m_cell_metrics->flx_pt_Jinv(cell_face_flx_pt) * inner_solution_gradient;


  compute_boundary_solution( inner_solution, inner_solution_gradient,
                             m_cell_metrics->flx_pt_coords(cell_face_flx_pt),
                             m_cell_metrics->flx_pt_unit_normal(cell_face_flx_pt)*m_cells_sf->flx_pt_sign(cell_face_flx_pt), 
                             boundary_solution, boundary_solution_gradient );
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                                   const ColVector_NDIM& coords,
                                                   const ColVector_NDIM& face_normal, 
                                                   RowVector_NEQS& boundary_solution )
{
  boundary_solution = inner_solution;
}


template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                                   const Matrix_NDIMxNEQS& inner_solution_gradient,
                                                   const ColVector_NDIM& coords,
                                                   const ColVector_NDIM& face_normal,
                                                   RowVector_NEQS& boundary_solution,
                                                   Matrix_NDIMxNEQS& boundary_solution_gradient )
{
  compute_boundary_solution(inner_solution, coords, face_normal, boundary_solution);
  boundary_solution_gradient = inner_solution_gradient;
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::set_face_solution( const Uint face_pt, 
                                           RowVector_NEQS& boundary_solution,
                                           Matrix_NDIMxNEQS& boundary_solution_gradient )
{
  const Uint pt = m_faces_space->connectivity()[m_face_idx][face_pt];
  cf3_assert(pt<bdry_solution()->size());
  cf3_assert(pt<bdry_solution_gradient()->size());
  cf3_assert(bdry_solution()->row_size() == NEQS);
  cf3_assert(bdry_solution_gradient()->row_size() == NEQS*NDIM);
  for (Uint eq=0; eq<NEQS; ++eq)
  {
    bdry_solution()->array()[pt][eq] = boundary_solution[eq];
    for( Uint d=0; d<NDIM; ++d)
    {
      bdry_solution_gradient()->array()[pt][eq+d*NEQS] = boundary_solution_gradient(d,eq);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3

#endif // cf3_dcm_core_BC_hpp
