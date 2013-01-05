// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_core_Metrics_hpp
#define cf3_sdm_core_Metrics_hpp

#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Connectivity.hpp"

#include "cf3/physics/MatrixTypes.hpp"

#include "cf3/sdm/core/LibCore.hpp"
#include "cf3/sdm/core/ShapeFunction.hpp"
#include "cf3/sdm/core/Reconstructions.hpp"

namespace cf3 {
namespace sdm {
namespace core {

/////////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM >
class sdm_core_API Metrics : public cf3::common::Component {

public: // types

  enum {NDIM = NB_DIM};

  typedef typename physics::MatrixTypes<NDIM>::ColVector_NDIM    ColVector_NDIM;
  typedef typename physics::MatrixTypes<NDIM>::Matrix_NDIMxNDIM  Matrix_NDIMxNDIM;

public: // functions

  /// Contructor
  /// @param name of the component
  Metrics ( const std::string& name );

  /// Virtual destructor
  virtual ~Metrics() {}

  /// Get the class name
  static std::string type_name () { return "Metrics"; }


  void setup_for_space(const Handle<mesh::Space const>& space);

  void compute_element(const Uint elem_idx);

  const mesh::ReconstructPoint&
    interpolation_from_cell_pts_to_flx_pt(const Uint flx_pt) const;
  
  const mesh::ReconstructPoint& 
    interpolation_from_sol_pts_to_flx_pt(const Uint flx_pt) const;

  const mesh::ReconstructPoint&
    set_zero(const Uint pt=0) const;

  const std::vector<mesh::ReconstructPoint>&
    zero_grad(const Uint pt=0) const;

  const mesh::ReconstructPoint&
    copy_pt(const Uint pt) const;

  const mesh::ReconstructPoint& 
    line_interpolation_from_flx_pts_to_sol_pt(const Uint sol_pt, const Uint line_direction) const;
  
  const mesh::ReconstructPoint&
    derivation_from_flx_pts_to_sol_pt(const Uint sol_pt, const Uint derivative_to) const;
  
  const mesh::ReconstructPoint&
    derivation_from_flx_pts_to_flx_pt(const Uint flx_pt, const Uint derivative_to) const;
  
  const mesh::ReconstructPoint&
    derivation_from_sol_pts_to_flx_pt(const Uint flx_pt, const Uint derivative_to) const;
      
  const std::vector< mesh::ReconstructPoint >&
    gradient_from_sol_pts_to_sol_pt(const Uint sol_pt) const;

  const std::vector< mesh::ReconstructPoint >&
    gradient_from_sol_pts_to_flx_pt(const Uint flx_pt) const;

  const std::vector< mesh::ReconstructPoint >&
    gradient_from_flx_pts_to_sol_pt(const Uint sol_pt) const;
  
  const ColVector_NDIM&
    sol_pt_coords(const Uint sol_pt) const;
  
//  const Matrix_NDIMxNDIM&
//    sol_pt_J(const Uint sol_pt) const;
  
//  const Matrix_NDIMxNDIM&
//    sol_pt_Jinv(const Uint sol_pt) const;
  
  const Real& 
    sol_pt_Jdet(const Uint sol_pt) const;
  
  const ColVector_NDIM& flx_pt_coords(const Uint flx_pt) const;
  
  const Matrix_NDIMxNDIM& 
    flx_pt_J(const Uint flx_pt) const;
  
  const Matrix_NDIMxNDIM& 
    flx_pt_Jinv(const Uint flx_pt) const;
  
  const Real& 
    flx_pt_Jdet(const Uint flx_pt) const;
  
  const ColVector_NDIM&
    flx_pt_unit_normal(const Uint flx_pt) const;
  
  const Real&
    flx_pt_Jvec_abs(const Uint flx_pt) const;
  
private:

  Uint m_elem_idx;

  Uint m_nb_faces;
  Uint m_nb_sol_pts;
  Uint m_nb_flx_pts;
  RealMatrix m_cell_coords;

  Handle<mesh::Cells const> m_cells;
  Handle<mesh::Space const> m_space;
  Handle<sdm::core::ShapeFunction const> m_sf;
  Handle<mesh::Field const> m_space_coords;

  std::vector< ColVector_NDIM >   m_sol_pt_coords;
  std::vector< Matrix_NDIMxNDIM > m_sol_pt_J;
  std::vector< Matrix_NDIMxNDIM > m_sol_pt_Jinv;
  std::vector< Real >             m_sol_pt_Jdet;
  
  std::vector< ColVector_NDIM >   m_flx_pt_coords;
  std::vector< Matrix_NDIMxNDIM > m_flx_pt_J;
  std::vector< Matrix_NDIMxNDIM > m_flx_pt_Jinv;
  std::vector< Real >             m_flx_pt_Jdet;
  std::vector< ColVector_NDIM >   m_flx_pt_unit_normal;
  std::vector< Real >             m_flx_pt_Jvec_abs;

  mesh::ReconstructPoint                                m_set_zero;
  std::vector< mesh::ReconstructPoint >                 m_zero_grad;
  std::vector< mesh::ReconstructPoint >                 m_copy_pt;
  std::vector< mesh::ReconstructPoint >                 m_interpolation_from_cell_pts_to_flx_pts;
  std::vector< mesh::ReconstructPoint >                 m_interpolation_from_sol_pts_to_flx_pts;
  std::vector< std::vector< mesh::ReconstructPoint > >  m_line_interpolation_from_flx_pts_to_sol_pts;
  std::vector< std::vector< mesh::ReconstructPoint > >  m_derivation_from_flx_pts_to_sol_pts;
  std::vector< std::vector< mesh::ReconstructPoint > >  m_derivation_from_flx_pts_to_flx_pts;
  std::vector< std::vector< mesh::ReconstructPoint > >  m_derivation_from_sol_pts_to_flx_pts;
  std::vector< std::vector< mesh::ReconstructPoint > >  m_derivation_from_sol_pts_to_sol_pts;

};

/////////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
Metrics<NB_DIM>::Metrics ( const std::string& name ) :
  cf3::common::Component(name)
{
}

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
void Metrics<NB_DIM>::setup_for_space(const Handle<mesh::Space const>& space)
{
  if (space != m_space)
  {
    m_space = space;
    cf3_assert(m_space);

    m_cells = space->support().handle<mesh::Cells>();
    cf3_assert(m_cells);

    m_sf = m_space->shape_function().handle<sdm::core::ShapeFunction>();
    cf3_assert(m_sf);

    m_space_coords = m_space->dict().coordinates().handle<mesh::Field>();
    cf3_assert(m_space_coords);

    m_nb_faces   = m_cells->element_type().nb_faces();
    m_nb_sol_pts = m_sf->nb_sol_pts();
    m_nb_flx_pts = m_sf->nb_flx_pts();

    m_cells->geometry_space().allocate_coordinates(m_cell_coords);

    m_sol_pt_coords .resize(m_nb_sol_pts);
    m_sol_pt_J      .resize(m_nb_sol_pts);
    m_sol_pt_Jinv   .resize(m_nb_sol_pts);
    m_sol_pt_Jdet   .resize(m_nb_sol_pts);

    m_flx_pt_coords      .resize(m_nb_flx_pts);
    m_flx_pt_J           .resize(m_nb_flx_pts);
    m_flx_pt_Jinv        .resize(m_nb_flx_pts);
    m_flx_pt_Jdet        .resize(m_nb_flx_pts);
    m_flx_pt_unit_normal .resize(m_nb_flx_pts);
    m_flx_pt_Jvec_abs    .resize(m_nb_flx_pts);
    
    m_interpolation_from_cell_pts_to_flx_pts .resize(m_nb_flx_pts);
    m_interpolation_from_sol_pts_to_flx_pts  .resize(m_nb_flx_pts);
    m_copy_pt                                .resize(m_nb_flx_pts);
    m_line_interpolation_from_flx_pts_to_sol_pts  .resize(m_nb_sol_pts, std::vector<mesh::ReconstructPoint>(NDIM));
    m_derivation_from_sol_pts_to_flx_pts     .resize(m_nb_flx_pts, std::vector<mesh::ReconstructPoint>(NDIM));
    m_derivation_from_flx_pts_to_sol_pts     .resize(m_nb_sol_pts, std::vector<mesh::ReconstructPoint >(NDIM));
    m_derivation_from_flx_pts_to_flx_pts     .resize(m_nb_flx_pts, std::vector<mesh::ReconstructPoint >(NDIM));
    m_derivation_from_sol_pts_to_sol_pts     .resize(m_nb_sol_pts, std::vector<mesh::ReconstructPoint >(NDIM));


    m_set_zero.m_N.resize(m_nb_flx_pts);
    m_set_zero.m_N.setZero();
    m_set_zero.m_pts.clear();

    m_zero_grad.resize(NDIM);
    for (Uint d=0; d<NDIM; ++d)
    {
      m_zero_grad[d].m_N.resize(m_nb_flx_pts);
      m_zero_grad[d].m_N.setZero();
      m_zero_grad[d].m_pts.clear();
    }

    for (Uint flx_pt=0; flx_pt<m_nb_flx_pts; ++flx_pt)
    {
      mesh::InterpolateInPoint::build_coefficients( m_interpolation_from_cell_pts_to_flx_pts[flx_pt],
                                                    m_sf->flx_pts().row(flx_pt),
                                                    m_cells->element_type().shape_function().handle<mesh::ShapeFunction>() );

      mesh::InterpolateInPoint::build_coefficients( m_interpolation_from_sol_pts_to_flx_pts[flx_pt],
                                                    m_sf->flx_pts().row(flx_pt),
                                                    m_sf->handle<mesh::ShapeFunction>() );

      for (Uint d=0; d<NDIM; ++d)
      {
        mesh::DerivativeInPoint::build_coefficients( m_derivation_from_sol_pts_to_flx_pts[flx_pt][d],
                                                     d ,
                                                     m_sf->flx_pts().row(flx_pt),
                                                     m_sf->handle<mesh::ShapeFunction>() );
        DerivativeInPointFromFlxPts::build_coefficients( m_derivation_from_flx_pts_to_flx_pts[flx_pt][d],
                                                         d,
                                                         m_sf->flx_pts().row(flx_pt),
                                                         m_sf );
      }

      m_copy_pt[flx_pt].m_N.resize(m_nb_flx_pts);
      m_copy_pt[flx_pt].m_N.setZero();
      m_copy_pt[flx_pt].m_N[flx_pt]=1.;
      m_copy_pt[flx_pt].m_pts.resize(1);
      m_copy_pt[flx_pt].m_pts[0]=flx_pt;
    }


    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      for (Uint d=0; d<NDIM; ++d)
      {
        InterpolateInPointFromFlxPts::build_coefficients( m_line_interpolation_from_flx_pts_to_sol_pts[sol_pt][d],
                                                          d,
                                                          m_sf->sol_pts().row(sol_pt),
                                                          m_sf );
        DerivativeInPointFromFlxPts::build_coefficients( m_derivation_from_flx_pts_to_sol_pts[sol_pt][d],
                                                         d,
                                                         m_sf->sol_pts().row(sol_pt),
                                                         m_sf );
        mesh::DerivativeInPoint::build_coefficients( m_derivation_from_sol_pts_to_sol_pts[sol_pt][d],
                                                     d,
                                                     m_sf->sol_pts().row(sol_pt),
                                                     m_sf->handle<mesh::ShapeFunction>() );
      }
    }
    m_elem_idx = math::Consts::uint_max(); // initialize it just to be different from zero
  }
}
////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
void Metrics<NB_DIM>::compute_element(const Uint elem_idx)
{
  if (elem_idx != m_elem_idx)
  {
    m_cells->geometry_space().put_coordinates(m_cell_coords,elem_idx);

    // Compute metrics in flux points
    for (Uint flx_pt=0; flx_pt<m_nb_flx_pts; ++flx_pt)
    {
      // Compute jacobian
      m_cells->element_type().compute_jacobian(m_sf->flx_pts().row(flx_pt), m_cell_coords,m_flx_pt_J[flx_pt]);
      m_flx_pt_Jinv[flx_pt] = m_flx_pt_J[flx_pt].inverse();
      m_flx_pt_Jdet[flx_pt] = m_flx_pt_J[flx_pt].determinant();
      cf3_always_assert(m_flx_pt_Jdet[flx_pt]>0);

      // Compute the unit-normal of the KSI, ETA or ZTA direction
      const Uint d = m_sf->flx_pt_dir(flx_pt);
      cf3_assert_desc(common::to_str(d)+"<"+common::to_str((Uint)NDIM),d<NDIM);
      m_flx_pt_unit_normal[flx_pt] = m_flx_pt_Jinv[flx_pt].col(d) * m_flx_pt_Jdet[flx_pt];
      m_flx_pt_Jvec_abs[flx_pt] = m_flx_pt_unit_normal[flx_pt].norm();
      cf3_always_assert(m_flx_pt_Jvec_abs[flx_pt]>0);
      m_flx_pt_unit_normal[flx_pt] /= m_flx_pt_Jvec_abs[flx_pt];

//      std::cout << "Jinv  = \n"<< m_flx_pt_Jinv[flx_pt]<<std::endl;
//      std::cout << "|Jvec["<<d<<"]| = " << m_flx_pt_Jvec_abs[flx_pt]<<std::endl;
//      std::cout << "normal = " << m_flx_pt_unit_normal[flx_pt] << std::endl;

      // Compute coordinates
      m_flx_pt_coords[flx_pt].setZero();
      boost_foreach (Uint n,  m_interpolation_from_cell_pts_to_flx_pts[flx_pt].used_points() )
      {
        m_flx_pt_coords[flx_pt] += m_interpolation_from_cell_pts_to_flx_pts[flx_pt].coeff(n) * m_cell_coords.row(n);
      }
    }

    // Compute metrics in solution
    mesh::Connectivity::ConstRow nodes = m_space->connectivity()[elem_idx];
    for (Uint sol_pt=0; sol_pt<m_nb_sol_pts; ++sol_pt)
    {
      // Compute jacobian
      // m_sol_pt_J[sol_pt] = m_cells->element_type().jacobian(m_sf->sol_pts().row(sol_pt), m_cell_coords);
      // m_sol_pt_Jinv[sol_pt] = m_sol_pt_J[sol_pt].inverse();
      // m_sol_pt_Jdet[sol_pt] = m_sol_pt_J[sol_pt].determinant();

      // Compute coordinates
      for (Uint d=0; d<NDIM; ++d)
        m_sol_pt_coords[sol_pt][d] = m_space_coords->array()[nodes[sol_pt]][d];

      m_sol_pt_Jdet[sol_pt] = m_cells->element_type().jacobian_determinant(m_sf->sol_pts().row(sol_pt), m_cell_coords);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::interpolation_from_cell_pts_to_flx_pt(const Uint flx_pt) const
{
  return m_interpolation_from_cell_pts_to_flx_pts[flx_pt];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint& 
  Metrics<NB_DIM>::interpolation_from_sol_pts_to_flx_pt(const Uint flx_pt) const
{
  return m_interpolation_from_sol_pts_to_flx_pts[flx_pt];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::copy_pt(const Uint pt) const
{
  return m_copy_pt[pt];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::set_zero(const Uint pt) const
{
  return m_set_zero;
}

template <Uint NB_DIM>
inline const std::vector<mesh::ReconstructPoint>&
  Metrics<NB_DIM>::zero_grad(const Uint pt) const
{
  return m_zero_grad;
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint& 
  Metrics<NB_DIM>::line_interpolation_from_flx_pts_to_sol_pt(const Uint sol_pt, const Uint line_direction) const
{
  return m_line_interpolation_from_flx_pts_to_sol_pts[sol_pt][line_direction];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::derivation_from_flx_pts_to_sol_pt(const Uint sol_pt, const Uint derivative_to) const
{
  return m_derivation_from_flx_pts_to_sol_pts[sol_pt][derivative_to];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::derivation_from_flx_pts_to_flx_pt(const Uint flx_pt, const Uint derivative_to) const
{
  return m_derivation_from_flx_pts_to_flx_pts[flx_pt][derivative_to];
}

template <Uint NB_DIM>
inline const mesh::ReconstructPoint&
  Metrics<NB_DIM>::derivation_from_sol_pts_to_flx_pt(const Uint flx_pt, const Uint derivative_to) const
{
  return m_derivation_from_sol_pts_to_flx_pts[flx_pt][derivative_to];
}

template <Uint NB_DIM>
inline const std::vector< mesh::ReconstructPoint >&
  Metrics<NB_DIM>::gradient_from_sol_pts_to_flx_pt(const Uint flx_pt) const
{
  return m_derivation_from_sol_pts_to_flx_pts[flx_pt];
}

template <Uint NB_DIM>
inline const std::vector< mesh::ReconstructPoint >&
  Metrics<NB_DIM>::gradient_from_sol_pts_to_sol_pt(const Uint sol_pt) const
{
  return m_derivation_from_sol_pts_to_sol_pts[sol_pt];
}

template <Uint NB_DIM>
inline const std::vector< mesh::ReconstructPoint >&
  Metrics<NB_DIM>::gradient_from_flx_pts_to_sol_pt(const Uint sol_pt) const
{
  return m_derivation_from_flx_pts_to_sol_pts[sol_pt];
}

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
inline const typename Metrics<NB_DIM>::ColVector_NDIM&
  Metrics<NB_DIM>::sol_pt_coords(const Uint sol_pt) const
{
  return m_sol_pt_coords[sol_pt];
}

//template <Uint NB_DIM>
//inline const typename Metrics<NB_DIM>::Matrix_NDIMxNDIM&
//  Metrics<NB_DIM>::sol_pt_J(const Uint sol_pt) const
//{
//  return m_sol_pt_J[sol_pt];
//}

//template <Uint NB_DIM>
//inline const typename Metrics<NB_DIM>::Matrix_NDIMxNDIM&
//  Metrics<NB_DIM>::sol_pt_Jinv(const Uint sol_pt) const
//{
//  return m_sol_pt_Jinv[sol_pt];
//}

template <Uint NB_DIM>
inline const Real&
  Metrics<NB_DIM>::sol_pt_Jdet(const Uint sol_pt) const
{
  return m_sol_pt_Jdet[sol_pt];
}

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
inline const typename Metrics<NB_DIM>::ColVector_NDIM&
  Metrics<NB_DIM>::flx_pt_coords(const Uint flx_pt) const
{
  return m_flx_pt_coords[flx_pt];
}

template <Uint NB_DIM>
inline const typename Metrics<NB_DIM>::Matrix_NDIMxNDIM&
  Metrics<NB_DIM>::flx_pt_J(const Uint flx_pt) const
{
  return m_flx_pt_J[flx_pt];
}

template <Uint NB_DIM>
inline const typename Metrics<NB_DIM>::Matrix_NDIMxNDIM&
  Metrics<NB_DIM>::flx_pt_Jinv(const Uint flx_pt) const
{
  return m_flx_pt_Jinv[flx_pt];
}

template <Uint NB_DIM>
inline const Real&
  Metrics<NB_DIM>::flx_pt_Jdet(const Uint flx_pt) const
{
  return m_flx_pt_Jdet[flx_pt];
}

template <Uint NB_DIM>
inline const typename Metrics<NB_DIM>::ColVector_NDIM&
  Metrics<NB_DIM>::flx_pt_unit_normal(const Uint flx_pt) const
{
  return m_flx_pt_unit_normal[flx_pt] ;
}

template <Uint NB_DIM>
inline const Real&
  Metrics<NB_DIM>::flx_pt_Jvec_abs(const Uint flx_pt) const
{
  return m_flx_pt_Jvec_abs[flx_pt];
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

#endif // cf3_sdm_core_Metrics_hpp
