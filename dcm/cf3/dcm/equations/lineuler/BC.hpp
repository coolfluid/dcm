// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the BCs of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_BC_hpp
#define cf3_dcm_equations_lineuler_BC_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"
#include "cf3/common/Option.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

/// @brief Base class for Boundary conditions that need the background flow
///        in Linearized Euler equations
/// @author Willem Deconinck
template <Uint NB_DIM, Uint NB_EQS>
class dcm_equations_lineuler_API BC : public dcm::core::BC<NB_DIM,NB_EQS>
{
public:
  BC(const std::string& name);
  virtual ~BC() {}
  static std::string type_name() { return "BC"; }

  virtual void compute_element(const Uint elem_idx);

  typedef typename dcm::core::BC<NB_DIM,NB_EQS>::ColVector_NDIM ColVector_NDIM;
  typedef typename dcm::core::BC<NB_DIM,NB_EQS>::RowVector_NEQS RowVector_NEQS;
protected:

  // These configuration variables are to compute the sound speed
  Real m_gamma;
  Real m_rho0;
  Real m_p0;
  ColVector_NDIM m_U0;
  Real m_c0;
  Handle< mesh::Field > m_background;

};

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM, Uint NB_EQS>
BC<NB_DIM,NB_EQS>::BC( const std::string& name ) :
  dcm::core::BC<NB_DIM,NB_EQS>(name),
  m_gamma(0.)
{
  this->options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  this->options().add("background",m_background).link_to(&m_background)
      .mark_basic()
      .description("background flow");
}

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM, Uint NB_EQS >
void BC<NB_DIM,NB_EQS>::compute_element( const Uint face_idx )
{
  this->m_face_idx = face_idx;
  // ---------------------------------------------------------------------------
  //                     COMPUTE CELL METRICS AND CONNECTIVITY
  // ---------------------------------------------------------------------------
  this->m_connected->compute_face(*this->m_faces,face_idx);
  cf3_assert(m_connected->is_bdry_face());
  this->m_cells = this->m_connected->cells()[LEFT].comp->template handle<mesh::Cells>();
  cf3_assert(this->m_cells);
  this->m_cells_space = this->m_fields->space(this->m_cells);
  this->m_cell_idx = this->m_connected->cells()[LEFT].idx;
  cf3_assert(m_cells_space);
  this->m_cells_sf = this->m_cells_space->shape_function().template handle<dcm::core::ShapeFunction>();
  cf3_assert(this->m_cells_sf);

  if  ( Handle<common::Component const> found = this->m_cells_space->get_child("metrics") )
  {
    this->m_metrics = const_cast<common::Component*>(found.get())->handle< dcm::core::Metrics<NB_DIM> >();
  }
  else
  {
    this->m_metrics = const_cast<mesh::Space*>(this->m_cells_space.get())->create_component< dcm::core::Metrics<NB_DIM> >("metrics");
    this->m_metrics->setup_for_space(this->m_cells_space);
  }
  this->m_cell_metrics = this->m_metrics->element(this->m_cell_idx);

  const std::vector<Uint>& cell_face_flx_pts =
    this->m_cells_sf->face_flx_pts( this->m_connected->cells_face_nb()[LEFT],
                                    this->m_connected->cells_orientation()[LEFT],
                                    this->m_connected->cells_rotation()[LEFT] );

  // Get flux-point characteristic solution
  std::vector< RowVector_NEQS > flx_pt_solution(this->m_cells_sf->nb_flx_pts());

  //CFinfo << "elem " << face_idx << "  interpolate sol" << CFendl;
  mesh::Connectivity::ConstRow nodes = this->m_cells_space->connectivity()[this->m_cell_idx];
  for (Uint f=0; f<this->m_cells_sf->nb_flx_pts(); ++f)
  {
    flx_pt_solution[f].setZero();
    boost_foreach( const Uint sol_pt, this->m_metrics->interpolation_from_sol_pts_to_flx_pt(f).used_points() )
    {
      const Real C = this->m_metrics->interpolation_from_sol_pts_to_flx_pt(f).coeff(sol_pt);
      for (Uint eq=0; eq<NB_EQS; ++eq)
        flx_pt_solution[f][eq] += C * this->solution()->array()[nodes[sol_pt]][eq];
    }
  }

  cf3_always_assert(m_background);
  // ---------------------------------------------------------------------------
  //                       COMPUTE BOUNDARY CONDITION
  // ---------------------------------------------------------------------------
  for (Uint face_pt=0; face_pt<this->m_nb_face_pts; ++face_pt)
  {
    const Uint cell_face_flx_pt = cell_face_flx_pts[face_pt];

    //CFinfo << "elem " << face_idx << "  interpolate bg" << CFendl;

    RowVector_NEQS background;
    background.setZero();
    boost_foreach( const Uint sol_pt, this->m_metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).used_points() )
    {
      const Real C = this->m_metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).coeff(sol_pt);
      for (Uint eq=0; eq<NB_EQS; ++eq)
      {
        background[eq] += C * this->m_background->array()[nodes[sol_pt]][eq];
      }
    }
    //CFinfo << "elem " << face_idx << "  extract bg" << CFendl;

    m_rho0 = background[0];
    for( int d=0; d<NB_DIM; ++d)
      m_U0[d] = background[d+1];
    m_p0 = background[3];
    m_c0 = std::sqrt(m_gamma*m_p0/m_rho0);

    //CFinfo << "elem " << face_idx << "  compute" << CFendl;
    this->compute_boundary_solution(cell_face_flx_pt, this->m_boundary_solution[face_pt], this->m_boundary_solution_gradient[face_pt]);

    this->set_face_solution(face_pt,this->m_boundary_solution[face_pt], this->m_boundary_solution_gradient[face_pt]);
  }
}

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_BC_hpp
