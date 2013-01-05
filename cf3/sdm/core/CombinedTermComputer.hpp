// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_core_CombinedTermComputer_hpp
#define cf3_sdm_core_CombinedTermComputer_hpp

#include "cf3/math/Consts.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/solver/TermComputer.hpp"
#include "cf3/physics/MatrixTypes.hpp"
#include "cf3/sdm/core/LibCore.hpp"
#include "cf3/sdm/core/CellConnectivity.hpp"
#include "cf3/sdm/core/ShapeFunction.hpp"
#include "cf3/sdm/core/Reconstructions.hpp"
#include "cf3/sdm/core/Metrics.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

template < typename TERM >
class sdm_core_API CombinedTermComputer : public solver::TermComputer {

public: // types

  enum {NDIM  = TERM::NDIM};
  enum {NEQS  = TERM::NEQS};
  enum {NVAR  = TERM::NVAR};
  enum {NGRAD = TERM::NGRAD};

  typedef physics::MatrixTypes<NDIM,NEQS,NVAR,NGRAD> Types;

  typedef typename TERM::DATA              PhysData;

public: // functions
  /// Contructor
  /// @param name of the component
  CombinedTermComputer ( const std::string& name );

  /// Virtual destructor
  virtual ~CombinedTermComputer() {}

  /// Get the class name
  static std::string type_name () { return TERM::type_name()+"Computer"; }

  virtual bool loop_cells(const Handle<mesh::Entities const>& cells);

  // virtual void compute_term(mesh::Field& term_field, mesh::Field& wave_speed) 
  // { 
  //   solver::TermComputer::compute_term(term_field,wave_speed); 
  // }

  virtual void compute_term(const Uint elem_idx, std::vector<RealVector>& term, std::vector<Real>& wave_speed);

private:

  Handle< mesh::Dictionary >  m_dict;
  Handle< TERM >              m_term;
  Handle< sdm::core::CellConnectivity >  m_connected;
  Handle< mesh::Cells const > m_cells;
  Handle< mesh::Space const > m_space;
  Handle< sdm::core::ShapeFunction const > m_sf;
  Handle< sdm::core::Metrics<NDIM> > m_metrics;

  Uint m_nb_faces;
  Uint m_nb_sol_pts;
  Uint m_nb_flx_pts;
  std::vector< typename Types::RowVector_NEQS > m_flx_pt_flux;
  std::vector<Real> m_flx_pt_wave_speed;
  std::vector< typename Types::ColVector_NDIM > m_sol_pt_wave_speed_vector;

  PhysData m_phys_data;
  PhysData m_neighbor_phys_data;

  std::vector< typename Types::RowVector_NEQS > m_sol_pt_flux_divergence;
  std::vector< typename Types::RowVector_NEQS > m_sol_pt_source;

  std::vector< typename Types::RowVector_NVAR > m_sol_pt_vars;
  std::vector< typename Types::RowVector_NVAR > m_flx_pt_vars;
  std::vector< typename Types::RowVector_NVAR > m_avg_flx_pt_vars;

  std::vector< typename Types::RowVector_NGRAD > m_sol_pt_gradvars;
  std::vector< typename Types::RowVector_NGRAD > m_flx_pt_gradvars;
  std::vector< typename Types::RowVector_NGRAD > m_avg_flx_pt_gradvars;

  std::vector< typename Types::Matrix_NDIMxNGRAD > m_sol_pt_gradvars_grad;
  std::vector< typename Types::Matrix_NDIMxNGRAD > m_flx_pt_gradvars_grad;
  std::vector< typename Types::Matrix_NDIMxNGRAD > m_avg_flx_pt_gradvars_grad;

  std::vector< typename Types::RowVector_NVAR > m_neighbor_flx_pt_vars;
  std::vector< typename Types::RowVector_NGRAD > m_neighbor_flx_pt_gradvars;
  std::vector< typename Types::Matrix_NDIMxNGRAD > m_neighbor_flx_pt_gradvars_grad;
  std::vector< typename Types::RowVector_NGRAD > m_flx_pt_gradvars_jump;
  std::vector< typename Types::RowVector_NGRAD > m_neighbor_flx_pt_gradvars_jump;
  std::vector< typename Types::Matrix_NDIMxNGRAD > m_LambdaL;
  std::vector< typename Types::Matrix_NDIMxNGRAD > m_LambdaR;

  Real m_convective_wave_speed;
  Real m_diffusive_wave_speed;
  typename Types::RowVector_NEQS m_convective_flux;
  typename Types::RowVector_NEQS m_diffusive_flux;


  Uint m_nb_face_pts;
  std::vector< std::vector<Uint> > m_face_pts;
  std::vector< std::vector<Uint> > m_neighbor_face_pts;
  std::vector< mesh::Space const* > m_neighbor_space;
  std::vector< Uint >               m_neighbor_elem_idx;

  Real m_alpha;
};

/////////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
CombinedTermComputer<TERM>::CombinedTermComputer ( const std::string& name ) :
  solver::TermComputer(name)
{
  options().add("term",m_term)
      .link_to(&m_term);

  options().add("alpha",-1.)
      .description("Damping coefficient in BR2 scheme for face-gradient computation\n"
                   "If negative, alpha = 1/order is used");
}

////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
bool CombinedTermComputer<TERM>::loop_cells(const Handle<mesh::Entities const>& cells)
{
  if ( is_null(cells->handle<mesh::Cells>()) )
    return false;

  if ( is_null(m_term) )
    throw common::SetupError(FromHere(), "term was not configured in "+uri().string() );

  m_dict = m_term->fields();
  cf3_assert(m_dict);

  m_cells = cells->handle<mesh::Cells>();
  cf3_assert(m_term);
  m_space = m_dict->space(m_cells);
  cf3_assert(m_space);
  m_sf = m_space->shape_function().handle<sdm::core::ShapeFunction>();
  cf3_assert(m_sf);

  // Set BR2 coefficient alpha to 1/order when alpha is negative
  m_alpha = options().template value<Real>("alpha");
  if (m_alpha < 0)
    m_alpha = 1./(Real)m_sf->order();


  // Create cellconnectivity
  if  ( Handle<Component const> found = m_space->get_child("cell_connectivity") )
  {
    m_connected = const_cast<Component*>(found.get())->handle< CellConnectivity >();
  }
  else
  {
    m_connected = const_cast<mesh::Space*>(m_space.get())->create_component< CellConnectivity >("cell_connectivity");
  }

  // Create and compute metrics for this space.
  // This means interpolation and derivation functions
  if  ( Handle<Component const> found = m_space->get_child("metrics") )
  {
    m_metrics = const_cast<Component*>(found.get())->handle< Metrics<NDIM> >();
  }
  else
  {
    m_metrics = const_cast<mesh::Space*>(m_space.get())->create_component< Metrics<NDIM> >("metrics");
    m_metrics->setup_for_space(m_space);
  }

  // Resize all vectors appropriately
  m_nb_sol_pts = m_sf->nb_sol_pts();
  m_nb_flx_pts = m_sf->nb_flx_pts();
  m_nb_faces = m_cells->element_type().nb_faces();
  m_nb_face_pts = m_sf->face_flx_pts(0,0,0).size();

  m_face_pts.resize( m_nb_faces, std::vector<Uint>(m_nb_face_pts) );
  m_neighbor_face_pts.resize( m_nb_faces, std::vector<Uint>(m_nb_face_pts) );
  m_neighbor_space.resize( m_nb_faces );
  m_neighbor_elem_idx.resize( m_nb_faces );

  m_sol_pt_vars.resize(m_nb_sol_pts);
  m_sol_pt_gradvars.resize(m_nb_sol_pts);
  m_sol_pt_gradvars_grad.resize(m_nb_sol_pts);
  m_sol_pt_flux_divergence.resize(m_nb_sol_pts,Types::RowVector_NEQS::Zero());
  m_sol_pt_wave_speed_vector.resize(m_nb_sol_pts);
  m_sol_pt_source.resize(m_nb_sol_pts,Types::RowVector_NEQS::Zero());

  m_flx_pt_flux.resize(m_nb_flx_pts);
  m_flx_pt_wave_speed.resize(m_nb_flx_pts);
  m_flx_pt_vars.resize(m_nb_flx_pts);
  m_flx_pt_gradvars.resize(m_nb_flx_pts);
  m_flx_pt_gradvars_grad.resize(m_nb_flx_pts);
  m_avg_flx_pt_vars.resize(m_nb_flx_pts);
  m_avg_flx_pt_gradvars.resize(m_nb_flx_pts);
  m_avg_flx_pt_gradvars_grad.resize(m_nb_flx_pts);
  m_LambdaL.resize(m_nb_flx_pts);
  m_LambdaR.resize(m_nb_flx_pts);
  m_flx_pt_gradvars_jump.resize(m_nb_flx_pts, Types::RowVector_NGRAD::Zero() );

  /// @warning It is assumed here that the neighbor element has the same element type as the current element
  ///          This means this won't work for neighbor cells with different order P
  m_neighbor_flx_pt_vars.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars_grad.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars_jump.resize(m_nb_flx_pts, Types::RowVector_NGRAD::Zero() );


  m_convective_flux.setZero();
  m_diffusive_flux.setZero();
  m_convective_wave_speed = 0.;
  m_diffusive_wave_speed  = 0.;
  return true;
}
////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
void CombinedTermComputer<TERM>::compute_term(const Uint elem_idx, std::vector<RealVector>& term, std::vector<Real>& wave_speed)
{

  // ---------------------------------------------------------------------------
  //                     COMPUTE CELL METRICS AND CONNECTIVITY
  // ---------------------------------------------------------------------------
  m_connected->compute_cell(*m_cells,elem_idx);
  m_metrics->compute_element(elem_idx);

  for (Uint face_nb=0; face_nb<m_nb_faces; ++face_nb)
  {
    m_neighbor_space[face_nb] = &m_dict->space(*m_connected->neighbor_cells()[face_nb].comp);
    m_neighbor_elem_idx[face_nb] = m_connected->neighbor_cells()[face_nb].idx;

    const std::vector<Uint>& face_flx_pts  = m_sf->face_flx_pts(
        face_nb,
        m_connected->orientations()[face_nb],
        m_connected->rotations()[face_nb] );

    if ( m_connected->is_bdry_face()[face_nb] )
    {
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        m_face_pts[face_nb][f]          = face_flx_pts[f];
        m_neighbor_face_pts[face_nb][f] = f;
      }
    }
    else
    {
      const std::vector<Uint>& neighbor_face_flx_pts = m_sf->face_flx_pts(
          m_connected->neighbour_face_nb()[face_nb],
          m_connected->neighbour_orientations()[face_nb],
          m_connected->neighbour_rotations()[face_nb] );
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        m_face_pts[face_nb][f]  = face_flx_pts[f];
        m_neighbor_face_pts[face_nb][f] = neighbor_face_flx_pts[f];
      }
    }
  }

  // ---------------------------------------------------------------------------
  //                           SET PHYSICAL CONSTANTS
  // ---------------------------------------------------------------------------
  m_term->set_phys_data_constants(m_phys_data);
  m_term->set_phys_data_constants(m_neighbor_phys_data);

  if (TERM::NVAR || TERM::NGRAD)
  {

    // ---------------------------------------------------------------------------
    //                   COMPUTE PHYSICAL VARIABLES AND GRADIENTS
    // ---------------------------------------------------------------------------
    for (Uint flx_pt=0; flx_pt<m_nb_flx_pts; ++flx_pt)
    {
      m_term->get_variables( /*in*/  *m_space, elem_idx,
                             /*in*/  m_metrics->flx_pt_coords(flx_pt),
                             /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(flx_pt),
                             /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(flx_pt),
                             /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
                             /*out*/ m_flx_pt_vars[flx_pt],
                             /*out*/ m_flx_pt_gradvars[flx_pt],
                             /*out*/ m_flx_pt_gradvars_grad[flx_pt] );
      m_avg_flx_pt_vars[flx_pt]          =  m_flx_pt_vars[flx_pt];
      m_avg_flx_pt_gradvars[flx_pt]      =  m_flx_pt_gradvars[flx_pt];
      m_avg_flx_pt_gradvars_grad[flx_pt] =  m_flx_pt_gradvars_grad[flx_pt];
    }

    // ---------------------------------------------------------------------------
    //                COMPUTE NEIGHBOR PHYSICAL VARIABLES AND GRADIENTS
    // ---------------------------------------------------------------------------
    for (Uint face_nb=0; face_nb<m_nb_faces; ++face_nb)
    {
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        Uint flx_pt = m_face_pts[face_nb][f];
        Uint neighbor_flx_pt = m_neighbor_face_pts[face_nb][f];

        if ( m_connected->is_bdry_face()[face_nb] )  // Gradient is copied from inside
        {
          m_term->get_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
                                 /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                 /*in*/  m_metrics->copy_pt(f),
                                 /*in*/  m_metrics->zero_grad(),
                                 /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
                                 /*out*/ m_neighbor_flx_pt_vars[flx_pt],
                                 /*out*/ m_neighbor_flx_pt_gradvars[flx_pt],
                                 /*out*/ m_neighbor_flx_pt_gradvars_grad[flx_pt] );

          m_neighbor_flx_pt_gradvars_grad[flx_pt] = m_flx_pt_gradvars_grad[flx_pt];
        }
        else // Gradient is computed in neighbor cell
        {
          m_term->get_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
                                 /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                 /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(neighbor_flx_pt),
                                 /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(neighbor_flx_pt),
                                 /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
                                 /*out*/ m_neighbor_flx_pt_vars[flx_pt],
                                 /*out*/ m_neighbor_flx_pt_gradvars[flx_pt],
                                 /*out*/ m_neighbor_flx_pt_gradvars_grad[flx_pt] );
        }

        m_avg_flx_pt_vars[flx_pt] += m_neighbor_flx_pt_vars[flx_pt];
        m_avg_flx_pt_vars[flx_pt] *= 0.5;
      }
    }
  }

  // ---------------------------------------------------------------------------
  //             COMPUTE FACE AVERAGE PHYSICAL VARIABLES AND GRADIENTS
  // ---------------------------------------------------------------------------
  //
  // Second approach of Bassi-Rebay (BR2).
  // -------------------------------------
  // Q_face = 1/2 * (Q_L + Q_R)
  // gradQ_face = 1/2 * (gradQ_L + gradQ_R) + alpha * 1/2 * (Lambda_L + Lambda_R)
  // LambdaL = grad ( dQL )  , computed in left cell
  // dQL = { (QR-QL)   in face
  //       {  0        elsewhere
  // LambdaR = grad ( dQR )  , computed in right cell
  // dQR = { (QL-QR)   in face
  //       {  0        elsewhere
  //
  // In a boundary face, Lambda_R = Lambda_L
  if (TERM::NGRAD)
  {
    for (Uint face_nb=0; face_nb<m_nb_faces; ++face_nb)
    {
      // Compute Left lifting operator ( Lambda_L )
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        m_flx_pt_gradvars_jump[m_face_pts[face_nb][f]] = m_neighbor_flx_pt_gradvars[m_face_pts[face_nb][f]] - m_flx_pt_gradvars[m_face_pts[face_nb][f]];
      }
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        typename Types::Matrix_NDIMxNGRAD& LambdaL = m_LambdaL[m_face_pts[face_nb][f]];
        LambdaL.setZero();
        for (Uint d=0; d<NDIM; ++d)
        {
          boost_foreach( Uint flx_pt, m_metrics->derivation_from_flx_pts_to_flx_pt(m_face_pts[face_nb][f],d).used_points() )
          {
            const Real C = m_metrics->derivation_from_flx_pts_to_flx_pt(m_face_pts[face_nb][f],d).coeff(flx_pt);
            for (Uint v=0; v<NGRAD; ++v)
              LambdaL(d,v) += m_flx_pt_gradvars_jump[flx_pt][v] * C;
          }
        }
        LambdaL = m_metrics->flx_pt_Jinv(m_face_pts[face_nb][f]) * LambdaL;
      }

      // Compute Right lifting operator ( Lambda_R )
      if ( m_connected->is_bdry_face()[face_nb] )
      {
        for (Uint f=0; f<m_nb_face_pts; ++f)
        {
          m_LambdaR[m_face_pts[face_nb][f]] = m_LambdaL[m_face_pts[face_nb][f]];
        }
      }
      else
      {
        for (Uint f=0; f<m_nb_face_pts; ++f)
        {
          m_neighbor_flx_pt_gradvars_jump[m_neighbor_face_pts[face_nb][f]] = - m_flx_pt_gradvars_jump[m_face_pts[face_nb][f]];
        }
        for (Uint f=0; f<m_nb_face_pts; ++f)
        {
          typename Types::Matrix_NDIMxNGRAD& LambdaR = m_LambdaR[m_face_pts[face_nb][f]];
          LambdaR.setZero();
          for (Uint d=0; d<NDIM; ++d)
          {
            boost_foreach( Uint flx_pt, m_metrics->derivation_from_flx_pts_to_flx_pt(m_neighbor_face_pts[face_nb][f],d).used_points() )
            {
              const Real C = m_metrics->derivation_from_flx_pts_to_flx_pt(m_neighbor_face_pts[face_nb][f],d).coeff(flx_pt);
              for (Uint v=0; v<NGRAD; ++v)
                LambdaR(d,v) += m_neighbor_flx_pt_gradvars_jump[flx_pt][v] * C;
            }
          }
          LambdaR = m_metrics->flx_pt_Jinv(m_face_pts[face_nb][f]) * LambdaR;
        }

        // set jumps back to zero in these flux points
        for (Uint f=0; f<m_nb_face_pts; ++f)
        {
          m_neighbor_flx_pt_gradvars_jump[m_neighbor_face_pts[face_nb][f]].setZero();
        }
      }
      // Average now
      for (Uint f=0; f<m_nb_face_pts; ++f) // loop over face flux points
      {
        m_avg_flx_pt_gradvars[m_face_pts[face_nb][f]] += m_neighbor_flx_pt_gradvars[m_face_pts[face_nb][f]];
        m_avg_flx_pt_gradvars[m_face_pts[face_nb][f]] *= 0.5;

        m_avg_flx_pt_gradvars_grad[m_face_pts[face_nb][f]] += m_neighbor_flx_pt_gradvars_grad[m_face_pts[face_nb][f]];
        m_avg_flx_pt_gradvars_grad[m_face_pts[face_nb][f]] += m_alpha*m_LambdaL[m_face_pts[face_nb][f]];
        m_avg_flx_pt_gradvars_grad[m_face_pts[face_nb][f]] += m_alpha*m_LambdaR[m_face_pts[face_nb][f]];
        m_avg_flx_pt_gradvars_grad[m_face_pts[face_nb][f]] *= 0.5;

        // set jumps to back to zero in these flux points
        m_flx_pt_gradvars_jump[m_face_pts[face_nb][f]].setZero();
      }
    }
  }

  if ( TERM::ENABLE_CONVECTION || TERM::ENABLE_DIFFUSION )
  {
    // -----------------------------------------------------------------------------------
    //   COMPUTE MAPPED FLUX IN INTERNAL FLUX POINTS IN DIRECTION 1_KSI, 1_ETA, OR 1_ZTA
    // -----------------------------------------------------------------------------------
    // computed wave_speed will be one of the following:
    //   ax*(dy*dz)  if  flux is in 1_KSI direction
    //   ay*(dx*dz)  if  flux is in 1_ETA direction
    //   az*(dx*dy)  if  flux is in 1_ZTA direction
    // ax , dx  are  wave speed  and  cell length in KSI direction
    // ay , dy  are  wave speed  and  cell length in ETA direction
    // az , dz  are  wave speed  and  cell length in ZTA direction
    boost_foreach (Uint flx_pt, m_sf->interior_flx_pts())
    {
      m_term->compute_phys_data( /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                 /*in*/  m_flx_pt_vars[flx_pt],
                                 /*in*/  m_flx_pt_gradvars[flx_pt],
                                 /*in*/  m_flx_pt_gradvars_grad[flx_pt],
                                 /*out*/ m_phys_data );

      // Compute fluxes projected on unit_normal
      m_term->compute_convective_flux( /*in*/  m_phys_data,
                                       /*in*/  m_metrics->flx_pt_unit_normal(flx_pt),
                                       /*out*/ m_convective_flux,
                                       /*out*/ m_convective_wave_speed );

      m_term->compute_diffusive_flux( /*in*/  m_phys_data,
                                      /*in*/  m_metrics->flx_pt_unit_normal(flx_pt),
                                      /*out*/ m_diffusive_flux,
                                      /*out*/ m_diffusive_wave_speed );

      m_flx_pt_flux[flx_pt]  = m_convective_flux;
      m_flx_pt_flux[flx_pt] -= m_diffusive_flux;

      m_flx_pt_wave_speed[flx_pt] = std::max(m_convective_wave_speed,m_diffusive_wave_speed);

      // Rescale flux and wavespeed because of transformation to mapped coordinate system
      m_flx_pt_flux[flx_pt]       *= m_metrics->flx_pt_Jvec_abs(flx_pt);
      m_flx_pt_wave_speed[flx_pt] *= m_metrics->flx_pt_Jvec_abs(flx_pt);
    }


    // --------------------------------------------------------------------------------
    //    COMPUTE MAPPED FLUX IN FACE FLUX POINTS, IN DIRECTION 1_KSI, 1_ETA, OR 1_ZTA
    // --------------------------------------------------------------------------------
    // computed wave_speed will be one of the following:
    //   ax*(dy*dz)  if  flux is in 1_KSI direction
    //   ay*(dx*dz)  if  flux is in 1_ETA direction
    //   az*(dx*dy)  if  flux is in 1_ZTA direction
    // ax , dx  are  wave speed  and  cell length in KSI direction
    // ay , dy  are  wave speed  and  cell length in ETA direction
    // az , dz  are  wave speed  and  cell length in ZTA direction
    for (Uint face_nb=0; face_nb<m_nb_faces; ++face_nb)
    {
      boost_foreach (Uint flx_pt, m_face_pts[face_nb])
      {
        if ( TERM::ENABLE_CONVECTION )
        {
          // compute physical data in flux points on this side
          m_term->compute_phys_data( /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_flx_pt_vars[flx_pt],
                                     /*in*/  m_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_phys_data );

          // compute physical data in flux points on other side
          m_term->compute_phys_data( /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_neighbor_flx_pt_vars[flx_pt],
                                     /*in*/  m_neighbor_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_neighbor_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_neighbor_phys_data );

          // Compute convective riemann-flux projected on OUTWARD unit_normal.
          // OUTWARD means that the sign of the unit-normal might have to be flipped, since the defined unit-normals
          // are always in positive KSI, ETA, ZTA direction.
          // The resulted flux must then also be added with a flipped sign
          m_term->compute_riemann_flux( /*in*/  m_phys_data,
                                        /*in*/  m_neighbor_phys_data,
                                        /*in*/  m_metrics->flx_pt_unit_normal(flx_pt)*m_sf->flx_pt_sign(flx_pt) ,
                                        /*out*/ m_convective_flux,
                                        /*out*/ m_convective_wave_speed );
          m_convective_flux *= m_sf->flx_pt_sign(flx_pt);
          m_flx_pt_flux[flx_pt] = m_convective_flux;
        }

        if ( TERM::ENABLE_DIFFUSION )
        {
          m_term->compute_phys_data( /*in*/  m_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_avg_flx_pt_vars[flx_pt],
                                     /*in*/  m_avg_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_avg_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_phys_data );

          m_term->compute_diffusive_flux( /*in*/  m_phys_data,
                                          /*in*/  m_metrics->flx_pt_unit_normal(flx_pt),
                                          /*out*/ m_diffusive_flux,
                                          /*out*/ m_diffusive_wave_speed );
          m_flx_pt_flux[flx_pt] -= m_diffusive_flux;
        }

        m_flx_pt_wave_speed[flx_pt] = std::max(m_convective_wave_speed,m_diffusive_wave_speed);

        // Rescale flux and wavespeed because of transformation to mapped coordinate system:
        m_flx_pt_flux[flx_pt]       *= m_metrics->flx_pt_Jvec_abs(flx_pt);
        m_flx_pt_wave_speed[flx_pt] *= m_metrics->flx_pt_Jvec_abs(flx_pt);
      }
    } // end for face_nb


    // ---------------------------------------------------------------------------
    //                  COMPUTE FLUX DIVERGENCE IN SOLUTION POINTS
    // ---------------------------------------------------------------------------
    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      // Computation of flux divergence in mapped coordinates
      m_sol_pt_flux_divergence[sol_pt].setZero();
      for (Uint d=0; d<NDIM; ++d)
      {
        boost_foreach( Uint flx_pt, m_metrics->derivation_from_flx_pts_to_sol_pt(sol_pt,d).used_points() )
        {

          const Real C = m_metrics->derivation_from_flx_pts_to_sol_pt(sol_pt,d).coeff(flx_pt);
          for (Uint eq=0; eq<NEQS; ++eq)
            m_sol_pt_flux_divergence[sol_pt][eq] += m_flx_pt_flux[flx_pt][eq] * C;
        }
      }
      // Transformation of flux divergence to physical coordinates
      m_sol_pt_flux_divergence[sol_pt] /= m_metrics->sol_pt_Jdet(sol_pt);
    }

    // ---------------------------------------------------------------------------
    //                   COMPUTE WAVE SPEEDS IN SOLUTION POINTS
    // ---------------------------------------------------------------------------
    // computed wave_speed will be of form:  ax/dx + ay/dy + az/dz
    // ax , dx  are  wave speed  and  cell length in KSI direction
    // ay , dy  are  wave speed  and  cell length in ETA direction
    // az , dz  are  wave speed  and  cell length in ZTA direction
    wave_speed.resize(m_nb_sol_pts);
    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      m_sol_pt_wave_speed_vector[sol_pt].setZero();
      wave_speed[sol_pt] = 0;
      // following assembles:  ax*(dy*dz) + ay*(dx*dz) + az*(dx*dy)
      for (Uint d=0; d<NDIM; ++d)
      {
        boost_foreach( Uint flx_pt, m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).used_points() )
        {
          m_sol_pt_wave_speed_vector[sol_pt][d]
              += m_flx_pt_wave_speed[flx_pt] * m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).coeff(flx_pt);
        }
        wave_speed[sol_pt] += m_sol_pt_wave_speed_vector[sol_pt][d];
      }
      // following divides:  ( ax*(dy*dz) + ay*(dx*dz) + az*(dx*dy)  )  /  (dx*dy*dz)
      wave_speed[sol_pt] /= m_metrics->sol_pt_Jdet(sol_pt);
    }
  } // end TERM::ENABLE_CONVECTION || TERM::ENABLE_DIFFUSION

  // ---------------------------------------------------------------------------
  //                        COMPUTE SOURCE IN SOLUTION POINTS
  // ---------------------------------------------------------------------------
  if (TERM::ENABLE_SOURCE)
  {
    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      if (TERM::NGRAD)
      {
        m_sol_pt_vars[sol_pt].setZero();
        m_sol_pt_gradvars[sol_pt].setZero();
        m_sol_pt_gradvars_grad[sol_pt].setZero();

        for (Uint d=0; d<NDIM; ++d)
        {
          boost_foreach( Uint flx_pt, m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).used_points() )
          {
            const Real C     = m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).coeff(flx_pt);
            const Real C_avg = C / NDIM;
            for (Uint v=0; v<NVAR; ++v)
            {
              m_sol_pt_vars[sol_pt][v] += m_avg_flx_pt_vars[flx_pt][v] * C_avg;
            }
            for (Uint v=0; v<NGRAD; ++v)
            {
              m_sol_pt_gradvars[sol_pt][v] += m_avg_flx_pt_gradvars[flx_pt][v] * C_avg;
              m_sol_pt_gradvars_grad[sol_pt](d,v) += m_avg_flx_pt_gradvars_grad[flx_pt](d,v) * C;
            }
          }
        }
      }
      m_term->compute_phys_data( /*in*/  m_metrics->sol_pt_coords(sol_pt),
                                 /*in*/  m_sol_pt_vars[sol_pt],
                                 /*in*/  m_sol_pt_gradvars[sol_pt],
                                 /*in*/  m_sol_pt_gradvars_grad[sol_pt],
                                 /*out*/ m_phys_data );
      m_term->compute_source( m_phys_data, m_sol_pt_source[sol_pt] );
    }
  }

  // ---------------------------------------------------------------------------
  //                           TERM = SRC + DIFF - CONV
  // ---------------------------------------------------------------------------
  term.resize(m_nb_sol_pts);
  for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
  {
    term[sol_pt].resize(NEQS);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      term[sol_pt][eq] = m_sol_pt_source[sol_pt][eq] - m_sol_pt_flux_divergence[sol_pt][eq];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

#endif // cf3_sdm_core_CombinedTermComputer_hpp
