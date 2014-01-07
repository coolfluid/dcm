// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_br2_BR2_hpp
#define cf3_sdm_br2_BR2_hpp

#include "cf3/math/Consts.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/solver/TermComputer.hpp"
#include "cf3/physics/MatrixTypes.hpp"
#include "cf3/sdm/br2/LibBR2.hpp"
#include "cf3/dcm/core/CellConnectivity.hpp"
#include "cf3/dcm/core/ShapeFunction.hpp"
#include "cf3/dcm/core/Reconstructions.hpp"
#include "cf3/dcm/core/Metrics.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace br2 {

////////////////////////////////////////////////////////////////////////////////

/// @brief Numerical computation of a term in a PDE
/// @note - It is assumed that all cells are of the same shape function
///       - flux points and solution points are aligned
template < typename TERM >
class sdm_br2_API BR2 : public solver::TermComputer {

public: // types

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  enum {NDIM  = TERM::NDIM};
  enum {NEQS  = TERM::NEQS};
  enum {NVAR  = TERM::NVAR};
  enum {NGRAD = TERM::NGRAD};

  typedef physics::MatrixTypes<NDIM,NEQS,NVAR,NGRAD> MatrixTypes;
  typedef typename MatrixTypes::RowVector_NEQS RowVector_NEQS;
  typedef typename MatrixTypes::RowVector_NVAR RowVector_NVAR;
  typedef typename MatrixTypes::ColVector_NDIM ColVector_NDIM;
  typedef typename MatrixTypes::RowVector_NGRAD RowVector_NGRAD;
  typedef typename MatrixTypes::Matrix_NDIMxNGRAD Matrix_NDIMxNGRAD;
  typedef typename TERM::DATA              PhysData;

public: // functions
  /// Contructor
  /// @param name of the component
  BR2 ( const std::string& name );

  /// Virtual destructor
  virtual ~BR2() {}

  /// Get the class name
  static std::string type_name () { return TERM::type_name()+"Computer"; }

  virtual bool loop_cells(const Handle<mesh::Entities const>& cells);

  virtual void compute_term(const Uint elem_idx, std::vector<RealVector>& term, std::vector<Real>& wave_speed);

private:

  Handle< mesh::Dictionary >  m_dict;
  Handle< mesh::Dictionary >  m_bdry_dict;
  Handle< TERM >              m_term;
  Handle< dcm::core::CellConnectivity >  m_connected;
  Handle< mesh::Cells const > m_cells;
  Handle< mesh::Space const > m_space;
  Handle< dcm::core::ShapeFunction const > m_sf;
  Handle< dcm::core::Metrics<NDIM> > m_metrics;
  Handle< dcm::core::ElementMetrics<NDIM> > m_element_metrics;
  std::vector< Handle< dcm::core::ElementMetrics<NDIM> > > m_neighbor_element_metrics;

  Uint m_nb_faces;
  Uint m_nb_sol_pts;
  Uint m_nb_flx_pts;
  std::vector< RowVector_NEQS > m_flx_pt_flux;
  std::vector<Real> m_flx_pt_wave_speed;
  std::vector< ColVector_NDIM > m_sol_pt_wave_speed_vector;

  PhysData m_phys_data;
  PhysData m_neighbor_phys_data;

  std::vector< RowVector_NEQS > m_sol_pt_flux_divergence;
  std::vector< RowVector_NEQS > m_sol_pt_source;

  std::vector< RowVector_NVAR > m_sol_pt_vars;
  std::vector< RowVector_NVAR > m_flx_pt_vars;
  std::vector< RowVector_NVAR > m_avg_flx_pt_vars;

  std::vector< RowVector_NGRAD > m_sol_pt_gradvars;
  std::vector< RowVector_NGRAD > m_flx_pt_gradvars;
  std::vector< RowVector_NGRAD > m_avg_flx_pt_gradvars;

  std::vector< Matrix_NDIMxNGRAD > m_sol_pt_gradvars_grad;
  std::vector< Matrix_NDIMxNGRAD > m_flx_pt_gradvars_grad;
  std::vector< Matrix_NDIMxNGRAD > m_avg_flx_pt_gradvars_grad;

  std::vector< RowVector_NVAR > m_neighbor_flx_pt_vars;
  std::vector< RowVector_NGRAD > m_neighbor_flx_pt_gradvars;
  std::vector< Matrix_NDIMxNGRAD > m_neighbor_flx_pt_gradvars_grad;
  std::vector< RowVector_NGRAD > m_flx_pt_gradvars_jump;
  std::vector< RowVector_NGRAD > m_neighbor_flx_pt_gradvars_jump;
  std::vector< Matrix_NDIMxNGRAD > m_LambdaL;
  std::vector< Matrix_NDIMxNGRAD > m_LambdaR;

  Real m_convective_wave_speed;
  Real m_diffusive_wave_speed;
  RowVector_NEQS m_convective_flux;
  RowVector_NEQS m_diffusive_flux;


  Uint m_nb_face_pts;
  std::vector< std::vector<Uint> > m_face_pts;
  std::vector< std::vector<Uint> > m_neighbor_face_pts;
  std::vector< mesh::Space const* > m_neighbor_space;
  std::vector< Uint >               m_neighbor_elem_idx;

  Real m_alpha;
};

/////////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
BR2<TERM>::BR2 ( const std::string& name ) :
  solver::TermComputer(name)
{
  m_term = create_static_component<TERM>("term");
  m_term->mark_basic();

  options().add("alpha",-1.)
      .description("Damping coefficient in BR2 scheme for face-gradient computation\n"
                   "If negative, alpha = 1/(P+1) is used");

  options().add("BR2",true)
      .description("Turn off neighbour cell influence");
}

////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
bool BR2<TERM>::loop_cells(const Handle<mesh::Entities const>& cells)
{
  if ( is_null(cells->handle<mesh::Cells>()) )
    return false;

  if ( is_null(m_term) )
    throw common::SetupError(FromHere(), "term was not configured in "+uri().string() );

  m_dict = m_term->fields();
  cf3_assert(m_dict);

  m_bdry_dict = m_term->bdry_fields();
  cf3_assert(m_bdry_dict);

  m_cells = cells->handle<mesh::Cells>();
  cf3_assert(m_term);
  m_space = m_dict->space(m_cells);
  cf3_assert(m_space);
  m_sf = m_space->shape_function().handle<dcm::core::ShapeFunction>();
  cf3_assert(m_sf);

  // Set BR2 coefficient alpha to 1/order when alpha is negative
  m_alpha = options().template value<Real>("alpha");
  if (m_alpha < 0)
  {
    m_alpha = 1./((Real)m_sf->order()+1.); // P0 --> solution is order 1
                                           //    --> alpha = 1.
  }


  // Create cellconnectivity
  if  ( Handle<Component const> found = m_space->get_child("cell_connectivity") )
  {
    m_connected = const_cast<Component*>(found.get())->handle< dcm::core::CellConnectivity >();
  }
  else
  {
    m_connected = const_cast<mesh::Space*>(m_space.get())->create_component< dcm::core::CellConnectivity >("cell_connectivity");
  }

  // Create and compute metrics for this space.
  // This means interpolation and derivation functions
  if  ( Handle<Component const> found = m_space->get_child("metrics") )
  {
    m_metrics = const_cast<Component*>(found.get())->handle< dcm::core::Metrics<NDIM> >();
  }
  else
  {
    m_metrics = const_cast<mesh::Space*>(m_space.get())->create_component< dcm::core::Metrics<NDIM> >("metrics");
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
  m_neighbor_element_metrics.resize( m_nb_faces );

  m_sol_pt_vars.resize(m_nb_sol_pts);
  m_sol_pt_gradvars.resize(m_nb_sol_pts);
  m_sol_pt_gradvars_grad.resize(m_nb_sol_pts);
  m_sol_pt_flux_divergence.resize(m_nb_sol_pts,RowVector_NEQS::Zero());
  m_sol_pt_wave_speed_vector.resize(m_nb_sol_pts);
  m_sol_pt_source.resize(m_nb_sol_pts,RowVector_NEQS::Zero());

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
  m_flx_pt_gradvars_jump.resize(m_nb_flx_pts, RowVector_NGRAD::Zero() );

  /// @warning It is assumed here that the neighbor element has the same element type as the current element
  ///          This means this won't work for neighbor cells with different order P
  m_neighbor_flx_pt_vars.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars_grad.resize(m_nb_flx_pts);
  m_neighbor_flx_pt_gradvars_jump.resize(m_nb_flx_pts, RowVector_NGRAD::Zero() );


  m_convective_flux.setZero();
  m_diffusive_flux.setZero();
  m_convective_wave_speed = 0.;
  m_diffusive_wave_speed  = 0.;

  return true;
}
////////////////////////////////////////////////////////////////////////////////

template <typename TERM >
void BR2<TERM>::compute_term(const Uint elem_idx, std::vector<RealVector>& term, std::vector<Real>& wave_speed)
{
  // ---------------------------------------------------------------------------
  //                     COMPUTE CELL METRICS AND CONNECTIVITY
  // ---------------------------------------------------------------------------
  m_connected->compute_cell(*m_cells,elem_idx);
  m_element_metrics = m_metrics->element(elem_idx);

  for (Uint face_nb=0; face_nb<m_nb_faces; ++face_nb)
  {
    m_neighbor_elem_idx[face_nb] = m_connected->neighbor_cells()[face_nb].idx;

    const std::vector<Uint>& face_flx_pts  = m_sf->face_flx_pts(
        face_nb,
        m_connected->orientations()[face_nb],
        m_connected->rotations()[face_nb] );

    if ( m_connected->is_bdry_face()[face_nb] )
    {
      m_neighbor_space[face_nb] = &m_bdry_dict->space(*m_connected->neighbor_cells()[face_nb].comp);
      for (Uint f=0; f<m_nb_face_pts; ++f)
      {
        m_face_pts[face_nb][f]          = face_flx_pts[f];
        m_neighbor_face_pts[face_nb][f] = f;
      }
    }
    else
    {
      m_neighbor_space[face_nb] = &m_dict->space(*m_connected->neighbor_cells()[face_nb].comp);
      m_neighbor_element_metrics[face_nb] = m_metrics->element( m_neighbor_elem_idx[face_nb] );
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
                             /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                             /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(flx_pt),
                             /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(flx_pt),
                             /*in*/  m_element_metrics->flx_pt_J(flx_pt),
                             /*in*/  m_element_metrics->flx_pt_Jinv(flx_pt),
                             /*in*/  m_element_metrics->flx_pt_Jdet(flx_pt),
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

        if ( m_connected->is_bdry_face()[face_nb] )
        {
          // Boundary values and Gradients are retrieved from fields, set by a BC
          m_term->get_bdry_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
                                      /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                                      /*in*/  m_metrics->copy_pt(f),
                                      /*in*/  m_metrics->copy_grad(f),
                                      /*in*/  m_element_metrics->flx_pt_J(flx_pt),
                                      /*in*/  m_element_metrics->flx_pt_Jinv(flx_pt),
                                      /*in*/  m_element_metrics->flx_pt_Jdet(flx_pt),
                                      /*out*/ m_neighbor_flx_pt_vars[flx_pt],
                                      /*out*/ m_neighbor_flx_pt_gradvars[flx_pt],
                                      /*out*/ m_neighbor_flx_pt_gradvars_grad[flx_pt] );
        }
        else
        {
          // Gradient is computed in neighbor cell
          m_term->get_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
                                 /*in*/  m_neighbor_element_metrics[face_nb]->flx_pt_coords(neighbor_flx_pt),
                                 /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(neighbor_flx_pt),
                                 /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(neighbor_flx_pt),
                                 /*in*/  m_neighbor_element_metrics[face_nb]->flx_pt_J(neighbor_flx_pt),
                                 /*in*/  m_neighbor_element_metrics[face_nb]->flx_pt_Jinv(neighbor_flx_pt),
                                 /*in*/  m_neighbor_element_metrics[face_nb]->flx_pt_Jdet(neighbor_flx_pt),
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
      for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
      {
        Uint flx_pt = m_face_pts[face_nb][face_pt];
        m_flx_pt_gradvars_jump[flx_pt] = m_neighbor_flx_pt_gradvars[flx_pt] - m_flx_pt_gradvars[flx_pt];
      }
      for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
      {
        Uint flx_pt = m_face_pts[face_nb][face_pt];
        Matrix_NDIMxNGRAD& LambdaL = m_LambdaL[flx_pt];
        LambdaL.setZero();
        for (Uint d=0; d<NDIM; ++d)
        {
          boost_foreach( Uint f, m_metrics->derivation_from_flx_pts_to_flx_pt(flx_pt,d).used_points() )
          {
            const Real C = m_metrics->derivation_from_flx_pts_to_flx_pt(flx_pt,d).coeff(f);
            for (Uint v=0; v<NGRAD; ++v)
              LambdaL(d,v) += m_flx_pt_gradvars_jump[f][v] * C;
          }
        }
        LambdaL = m_element_metrics->flx_pt_Jinv(flx_pt) * LambdaL;
      }

      // Compute Right lifting operator ( Lambda_R )
      if ( m_connected->is_bdry_face()[face_nb] )
      {
        for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
        {
          Uint flx_pt = m_face_pts[face_nb][face_pt];
          m_LambdaR[flx_pt] = m_LambdaL[flx_pt];
        }
      }
      else
      {
        for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
        {
          Uint flx_pt = m_face_pts[face_nb][face_pt];
          Uint neighbour_flx_pt = m_neighbor_face_pts[face_nb][face_pt];
          m_neighbor_flx_pt_gradvars_jump[neighbour_flx_pt] = - m_flx_pt_gradvars_jump[flx_pt];
        }
        for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
        {
          Uint flx_pt = m_face_pts[face_nb][face_pt];
          Uint neighbour_flx_pt = m_neighbor_face_pts[face_nb][face_pt];
          Matrix_NDIMxNGRAD& LambdaR = m_LambdaR[flx_pt];
          LambdaR.setZero();
          for (Uint d=0; d<NDIM; ++d)
          {
            boost_foreach( Uint f, m_metrics->derivation_from_flx_pts_to_flx_pt(neighbour_flx_pt,d).used_points() )
            {
              const Real C = m_metrics->derivation_from_flx_pts_to_flx_pt(neighbour_flx_pt,d).coeff(f);
              for (Uint v=0; v<NGRAD; ++v)
                LambdaR(d,v) += m_neighbor_flx_pt_gradvars_jump[f][v] * C;
            }
          }
          LambdaR = m_neighbor_element_metrics[face_nb]->flx_pt_Jinv(neighbour_flx_pt) * LambdaR;
        }

        // set jumps back to zero in these flux points
        for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt)
        {
          Uint neighbour_flx_pt = m_neighbor_face_pts[face_nb][face_pt];
          m_neighbor_flx_pt_gradvars_jump[neighbour_flx_pt].setZero();
        }
      }
      // Average now
      for (Uint face_pt=0; face_pt<m_nb_face_pts; ++face_pt) // loop over face flux points
      {
        Uint flx_pt = m_face_pts[face_nb][face_pt];
        m_avg_flx_pt_gradvars[flx_pt] += m_neighbor_flx_pt_gradvars[flx_pt];
        m_avg_flx_pt_gradvars[flx_pt] *= 0.5;
        m_avg_flx_pt_gradvars_grad[flx_pt] += m_neighbor_flx_pt_gradvars_grad[flx_pt];
        m_avg_flx_pt_gradvars_grad[flx_pt] += m_alpha*m_LambdaL[flx_pt];
        m_avg_flx_pt_gradvars_grad[flx_pt] += m_alpha*m_LambdaR[flx_pt];
        m_avg_flx_pt_gradvars_grad[flx_pt] *= 0.5;

        // set jumps to back to zero in these flux points
        m_flx_pt_gradvars_jump[flx_pt].setZero();
      }
    }
  }

  if ( TERM::ENABLE_CONVECTION || TERM::ENABLE_DIFFUSION )
  {
    // -----------------------------------------------------------------------------------
    //   COMPUTE MAPPED FLUX IN INTERNAL FLUX POINTS IN DIRECTION 1_KSI, 1_ETA, OR 1_ZTA
    // -----------------------------------------------------------------------------------
    // computed wave_speed will be the maximum of convective and diffusive wave speeds
    //                    KSI           ETA         ZTA
    //   - convective:  ax / dx       ay / dy     az / dz
    //   - diffusive:   mu / dx^2     mu / dy^2   mu / dz^2
    boost_foreach (Uint flx_pt, m_sf->interior_flx_pts())
    {
      m_term->compute_phys_data( /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                                 /*in*/  m_flx_pt_vars[flx_pt],
                                 /*in*/  m_flx_pt_gradvars[flx_pt],
                                 /*in*/  m_flx_pt_gradvars_grad[flx_pt],
                                 /*out*/ m_phys_data );

      m_flx_pt_flux[flx_pt].setZero();
      m_convective_wave_speed = 0.;
      m_diffusive_wave_speed = 0.;

      if ( TERM::ENABLE_CONVECTION )
      {
        // Compute fluxes projected on unit_normal
        m_term->compute_convective_flux( /*in*/  m_phys_data,
                                         /*in*/  m_element_metrics->flx_pt_unit_normal(flx_pt),
                                         /*out*/ m_convective_flux,
                                         /*out*/ m_convective_wave_speed );
        m_flx_pt_flux[flx_pt] += m_convective_flux;

        // Convective wavespeed --> a/dx
        m_convective_wave_speed *= m_element_metrics->flx_pt_Jvec_abs(flx_pt);
        m_convective_wave_speed /= m_element_metrics->flx_pt_Jdet(flx_pt);
      }

      if ( TERM::ENABLE_DIFFUSION )
      {
        // Compute fluxes projected on unit_normal
        m_term->compute_diffusive_flux( /*in*/  m_phys_data,
                                        /*in*/  m_element_metrics->flx_pt_unit_normal(flx_pt),
                                        /*out*/ m_diffusive_flux,
                                        /*out*/ m_diffusive_wave_speed );
        m_flx_pt_flux[flx_pt] -= m_diffusive_flux;

        // Diffusive wavespeed  --> mu/dx2
        m_diffusive_wave_speed *= m_element_metrics->flx_pt_Jvec_abs(flx_pt)*m_element_metrics->flx_pt_Jvec_abs(flx_pt);
        m_diffusive_wave_speed /= m_element_metrics->flx_pt_Jdet(flx_pt)*m_element_metrics->flx_pt_Jdet(flx_pt);
      }

      // Rescale flux because of transformation to mapped coordinate system
      m_flx_pt_flux[flx_pt]       *= m_element_metrics->flx_pt_Jvec_abs(flx_pt);

      // Take max of both wave speeds
      m_flx_pt_wave_speed[flx_pt] = std::max(m_convective_wave_speed,m_diffusive_wave_speed);
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
        m_flx_pt_flux[flx_pt].setZero();
        if ( TERM::ENABLE_CONVECTION )
        {
          // compute physical data in flux points on this side
          m_term->compute_phys_data( /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_flx_pt_vars[flx_pt],
                                     /*in*/  m_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_phys_data );

          // compute physical data in flux points on other side
          m_term->compute_phys_data( /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_neighbor_flx_pt_vars[flx_pt],
                                     /*in*/  m_neighbor_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_neighbor_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_neighbor_phys_data );

          // Compute convective riemann-flux projected on OUTWARD unit_normal.
          // OUTWARD means that the sign of the unit-normal MIGHT have to be flipped, since the defined unit-normals
          // are always in positive KSI, ETA, ZTA direction.
          // The resulted flux must then also be added with a flipped sign
          m_term->compute_riemann_flux( /*in*/  m_phys_data,
                                        /*in*/  m_neighbor_phys_data,
                                        /*in*/  m_element_metrics->flx_pt_unit_normal(flx_pt)*m_sf->flx_pt_sign(flx_pt) ,
                                        /*out*/ m_convective_flux,
                                        /*out*/ m_convective_wave_speed );
          m_convective_flux *= m_sf->flx_pt_sign(flx_pt);
          m_flx_pt_flux[flx_pt] += m_convective_flux;
        }

        if ( TERM::ENABLE_DIFFUSION )
        {
          m_term->compute_phys_data( /*in*/  m_element_metrics->flx_pt_coords(flx_pt),
                                     /*in*/  m_avg_flx_pt_vars[flx_pt],
                                     /*in*/  m_avg_flx_pt_gradvars[flx_pt],
                                     /*in*/  m_avg_flx_pt_gradvars_grad[flx_pt],
                                     /*out*/ m_phys_data );

          m_term->compute_diffusive_flux( /*in*/  m_phys_data,
                                          /*in*/  m_element_metrics->flx_pt_unit_normal(flx_pt)*m_sf->flx_pt_sign(flx_pt),
                                          /*out*/ m_diffusive_flux,
                                          /*out*/ m_diffusive_wave_speed );
          m_diffusive_flux *= m_sf->flx_pt_sign(flx_pt);
          m_flx_pt_flux[flx_pt] -= m_diffusive_flux;
        }

        // Rescale flux because of transformation to mapped coordinate system
        m_flx_pt_flux[flx_pt] *= m_element_metrics->flx_pt_Jvec_abs(flx_pt);
        // Convective wavespeed --> a/dx
        m_convective_wave_speed *= m_element_metrics->flx_pt_Jvec_abs(flx_pt);
        m_convective_wave_speed /= m_element_metrics->flx_pt_Jdet(flx_pt);

        // Diffusive wavespeed  --> mu/dx2
        m_diffusive_wave_speed *= m_element_metrics->flx_pt_Jvec_abs(flx_pt)*m_element_metrics->flx_pt_Jvec_abs(flx_pt);
        m_diffusive_wave_speed /= m_element_metrics->flx_pt_Jdet(flx_pt)*m_element_metrics->flx_pt_Jdet(flx_pt);

        // Take max of both
        m_flx_pt_wave_speed[flx_pt] = std::max(m_convective_wave_speed,m_diffusive_wave_speed);
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
      m_sol_pt_flux_divergence[sol_pt] /= m_element_metrics->sol_pt_Jdet(sol_pt);
    }

    // ---------------------------------------------------------------------------
    //                   COMPUTE WAVE SPEEDS IN SOLUTION POINTS
    // ---------------------------------------------------------------------------
    wave_speed.resize(m_nb_sol_pts);
    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      m_sol_pt_wave_speed_vector[sol_pt].setZero();
      wave_speed[sol_pt] = 0;
      for (Uint d=0; d<NDIM; ++d)
      {
        boost_foreach( Uint flx_pt, m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).used_points() )
        {
          m_sol_pt_wave_speed_vector[sol_pt][d]
              += m_flx_pt_wave_speed[flx_pt] * m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).coeff(flx_pt);
        }
        wave_speed[sol_pt] += m_sol_pt_wave_speed_vector[sol_pt][d];
      }
    }
  } // end TERM::ENABLE_CONVECTION || TERM::ENABLE_DIFFUSION

  // ---------------------------------------------------------------------------
  //                        COMPUTE SOURCE IN SOLUTION POINTS
  // ---------------------------------------------------------------------------
  if (TERM::ENABLE_SOURCE)
  {
    if (TERM::NGRAD)
    {
      // Transform gradient in flux points to KSI,ETA,ZTA coordinates
      for (Uint flx_pt=0; flx_pt<m_nb_flx_pts; ++flx_pt)
      {
        m_avg_flx_pt_gradvars_grad[flx_pt] = m_element_metrics->flx_pt_J(flx_pt) * m_avg_flx_pt_gradvars_grad[flx_pt];
      }
    }

    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      m_sol_pt_vars[sol_pt].setZero();
      if (TERM::NGRAD)
      {
        m_sol_pt_gradvars[sol_pt].setZero();
        m_sol_pt_gradvars_grad[sol_pt].setZero();
      }
      for (Uint d=0; d<NDIM; ++d)
      {
        boost_foreach( Uint flx_pt, m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).used_points() )
        {
          const Real C     = m_metrics->line_interpolation_from_flx_pts_to_sol_pt(sol_pt,d).coeff(flx_pt);
          const Real C_avg = C / static_cast<Real>(NDIM);
          for (Uint v=0; v<NVAR; ++v)
          {
            m_sol_pt_vars[sol_pt][v] += m_avg_flx_pt_vars[flx_pt][v] * C_avg;
          }
          for (Uint v=0; v<NGRAD; ++v)
          {
            m_sol_pt_gradvars[sol_pt][v]        += C_avg * m_avg_flx_pt_gradvars[flx_pt][v];
            m_sol_pt_gradvars_grad[sol_pt](d,v) += C     * m_avg_flx_pt_gradvars_grad[flx_pt](d,v);
          }
        }
      }
      // Transform interpolated gradient in solution points to X,Y,Z coordinates
      m_sol_pt_gradvars_grad[sol_pt] = m_element_metrics->sol_pt_Jinv(sol_pt) * m_sol_pt_gradvars_grad[sol_pt];

      m_term->compute_phys_data( /*in*/  m_element_metrics->sol_pt_coords(sol_pt),
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

} // br2
} // sdm
} // cf3

#endif // cf3_sdm_br2_BR2_hpp
