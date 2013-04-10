// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_tools_ComputeFieldGradientBR2_hpp
#define cf3_sdm_tools_ComputeFieldGradientBR2_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/mesh/MeshTransformer.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/sdm/tools/LibTools.hpp"
#include "cf3/sdm/core/Metrics.hpp"
#include "cf3/sdm/core/CellConnectivity.hpp"


////////////////////////////////////////////////////////////////////////////////

//Forward declarations
namespace cf3 {
  namespace mesh {
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace tools {

//////////////////////////////////////////////////////////////////////////////

/// @brief Compute the gradient of a field
/// @author Willem Deconinck
///
/// If field contains variables (in order) "a , b , c",
/// then field_gradient will contain variables (in order, e.g. 2D)
/// "da/dx  db/dx  dc/dx  da/dy db/dy dc/dy"
///
/// Shapefunctions are used to compute the gradient.
///
/// @note In the case of a P0 shape-function,
///       the gradient will be thus wrongly calculated to
///       be zero, as variables are piece-wise constant.
///
class sdm_tools_API ComputeFieldGradientBR2 : public mesh::MeshTransformer
{   
public: // functions
  
  /// constructor
  ComputeFieldGradientBR2( const std::string& name );
  
  /// Gets the Class name
  static std::string type_name() { return "ComputeFieldGradientBR2"; }

  virtual void execute();
  
  template< Uint NDIM >
  void compute_gradient();

private: // data

  Handle<mesh::Field const> m_field;
  Handle<mesh::Field>       m_field_gradient;
  std::vector<Real>   m_normal;
  Real m_alpha;


}; // end ComputeFieldGradientBR2

////////////////////////////////////////////////////////////////////////////////

template< Uint NDIM >
void ComputeFieldGradientBR2::compute_gradient()
{
  using namespace cf3::sdm::core;
  using namespace cf3::mesh;

  Handle< Metrics<NDIM> > metrics;
  Handle< CellConnectivity >  connected;

  const Uint nb_vars = m_field->row_size();

  // Loop over patches
  boost_foreach(const Handle<Space>& space, m_field->spaces())
  {
    if (space->shape_function().dimensionality() == NDIM) // if volume element
    {
      Handle<Cells> cells = space->support().handle<Cells>();
      Handle<const sdm::core::ShapeFunction> sf = space->shape_function().handle<sdm::core::ShapeFunction>();

      // Set BR2 coefficient alpha to 1/order when alpha is negative
      m_alpha = options().template value<Real>("alpha");
      if (m_alpha < 0)
      {
        m_alpha = 1./((Real)sf->order()+1.); // P0 --> solution is order 1
                                             //    --> alpha = 1.
      }


      if  ( Handle<Component const> found = space->get_child("metrics") )
      {
        metrics = const_cast<Component*>(found.get())->handle< Metrics<NDIM> >();
      }
      else
      {
        metrics = space->create_component< Metrics<NDIM> >("metrics");
        metrics->setup_for_space(space);
      }
      if  ( Handle<Component const> found = space->get_child("cell_connectivity") )
      {
        connected = const_cast<Component*>(found.get())->handle< CellConnectivity >();
      }
      else
      {
        connected = const_cast<mesh::Space*>(space.get())->create_component< CellConnectivity >("cell_connectivity");
      }

      const Uint nb_flx_pts = sf->nb_flx_pts();
      const Uint nb_sol_pts = sf->nb_sol_pts();
      const Uint nb_faces = cells->element_type().nb_faces();
      const Uint nb_face_pts = sf->face_flx_pts(0,0,0).size();

      std::vector< Space const* > neighbor_space(nb_faces);
      std::vector< Uint > neighbor_elem_idx(nb_faces);
      std::vector< std::vector<Uint> > face_pts(nb_faces, std::vector<Uint>(nb_face_pts));
      std::vector< std::vector<Uint> > neighbor_face_pts(nb_faces, std::vector<Uint>(nb_face_pts));

      std::vector<RealVector> flx_pt_vars(nb_flx_pts,RealVector(nb_vars));
      std::vector<RealVector> neighbor_flx_pt_vars(nb_flx_pts,RealVector(nb_vars));
      std::vector<RealVector> avg_flx_pt_vars(nb_flx_pts,RealVector(nb_vars));
      std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > flx_pt_grad(nb_flx_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
      std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > neighbor_flx_pt_grad(nb_flx_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
      std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > avg_flx_pt_grad(nb_flx_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
      std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > LambdaL(nb_flx_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
      std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > LambdaR(nb_flx_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
      std::vector<RealVector> jump(nb_flx_pts, RealVector(nb_vars));
      std::vector<RealVector> neighbor_jump(nb_flx_pts, RealVector(nb_vars));


      Handle<Space const> grad_space = m_field_gradient->space(*cells).handle<Space>();
      Handle<const mesh::ShapeFunction> grad_sf = grad_space->shape_function().handle<mesh::ShapeFunction>();
      const Uint nb_grad_pts = grad_sf->nb_nodes();

      std::vector< std::vector< mesh::ReconstructPoint > >  line_interpolation_from_flx_pts_to_grad_pt;
      line_interpolation_from_flx_pts_to_grad_pt.resize(nb_grad_pts, std::vector<mesh::ReconstructPoint>(NDIM));
      for (Uint grad_pt=0; grad_pt<nb_grad_pts; ++grad_pt)
      {
        for (Uint d=0; d<NDIM; ++d)
        {
          InterpolateInPointFromFlxPts::build_coefficients( line_interpolation_from_flx_pts_to_grad_pt[grad_pt][d],
                                                            d,
                                                            grad_sf->local_coordinates().row(grad_pt),
                                                            sf );
        }
      }


//      std::vector< std::vector< mesh::ReconstructPoint > >  derivation_from_flx_pts_to_grad_pt;
//      derivation_from_flx_pts_to_grad_pt.resize(nb_grad_pts, std::vector<mesh::ReconstructPoint >(NDIM));
//      for (Uint grad_pt=0; grad_pt<nb_grad_pts; ++grad_pt)
//      {
//        for (Uint d=0; d<NDIM; ++d)
//        {
//          DerivativeInPointFromFlxPts::build_coefficients( interpolation_from_flx_pts_to_grad_pt[grad_pt][d],
//                                                           d,
//                                                           grad_sf->local_coordinates().row(grad_pt),
//                                                           sf );
//        }
//      }


      // Loop over cells
      Uint nb_elems = space->size();
      for (Uint e=0; e<nb_elems; ++e)
      {
        // Compute cell metrics
        connected->compute_cell(*cells,e);
        metrics->compute_element(e);

        for (Uint face_nb=0; face_nb<nb_faces; ++face_nb)
        {
          neighbor_space[face_nb] = &m_field->dict().space(*connected->neighbor_cells()[face_nb].comp);
          neighbor_elem_idx[face_nb] = connected->neighbor_cells()[face_nb].idx;

          const std::vector<Uint>& face_flx_pts  = sf->face_flx_pts(
              face_nb,
              connected->orientations()[face_nb],
              connected->rotations()[face_nb] );

          if ( connected->is_bdry_face()[face_nb] )
          {
            for (Uint f=0; f<nb_face_pts; ++f)
            {
              face_pts[face_nb][f]          = face_flx_pts[f];
              neighbor_face_pts[face_nb][f] = f;
            }
          }
          else
          {
            const std::vector<Uint>& neighbor_face_flx_pts = sf->face_flx_pts(
                connected->neighbour_face_nb()[face_nb],
                connected->neighbour_orientations()[face_nb],
                connected->neighbour_rotations()[face_nb] );
            for (Uint f=0; f<nb_face_pts; ++f)
            {
              face_pts[face_nb][f]  = face_flx_pts[f];
              neighbor_face_pts[face_nb][f] = neighbor_face_flx_pts[f];
            }
          }

        }

        // Interpolate field to flux points
        for (Uint flx_pt=0; flx_pt<nb_flx_pts; ++flx_pt)
        {
          mesh::Connectivity::ConstRow nodes = space->connectivity()[e];

          const mesh::ReconstructPoint& interpolation = metrics->interpolation_from_sol_pts_to_flx_pt(flx_pt);
          flx_pt_vars[flx_pt].setZero();
          boost_foreach( Uint n, interpolation.used_points() )
          {
            const Real C = interpolation.coeff(n);
            for (Uint v=0; v<nb_vars; ++v)
            {
              flx_pt_vars[flx_pt][v] += C * m_field->array()[nodes[n]][v];
            }
          }
//          CFinfo << "var = " << flx_pt_vars[flx_pt] << CFendl;
          avg_flx_pt_vars[flx_pt] = flx_pt_vars[flx_pt];

          const std::vector<mesh::ReconstructPoint>& gradient = metrics->gradient_from_sol_pts_to_flx_pt(flx_pt);
          flx_pt_grad[flx_pt].setZero();
          for (Uint d=0; d<NDIM; ++d)
          {
            boost_foreach( Uint n, gradient[d].used_points() )
            {
              const Real C = gradient[d].coeff(n);
              for (Uint v=0; v<nb_vars; ++v)
              {
                flx_pt_grad[flx_pt](d,v) += C * m_field->array()[nodes[n]][v];
              }
            }
          }
//          flx_pt_grad[flx_pt] = metrics->flx_pt_Jinv(flx_pt)*flx_pt_grad[flx_pt];
//          CFinfo << flx_pt_grad[flx_pt].transpose() << CFendl;
          avg_flx_pt_grad[flx_pt] = flx_pt_grad[flx_pt];

//          CFinfo << flx_pt_vars[flx_pt] << CFendl;
          //          m_term->get_variables( /*in*/  *m_space, elem_idx,
          //                                 /*in*/  m_metrics->flx_pt_coords(flx_pt),
          //                                 /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(flx_pt),
          //                                 /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(flx_pt),
          //                                 /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
          //                                 /*out*/ m_flx_pt_vars[flx_pt],
          //                                 /*out*/ m_flx_pt_gradvars[flx_pt],
          //                                 /*out*/ m_flx_pt_gradvars_grad[flx_pt] );

        }

        for (Uint face_nb=0; face_nb<nb_faces; ++face_nb)
        {
          for (Uint f=0; f<nb_face_pts; ++f)
          {
            Uint flx_pt = face_pts[face_nb][f];
            Uint neighbor_flx_pt = neighbor_face_pts[face_nb][f];

            if ( connected->is_bdry_face()[face_nb] )  // Gradient is copied from inside
            {
              mesh::Connectivity::ConstRow nodes = neighbor_space[face_nb]->connectivity()[neighbor_elem_idx[face_nb]];

              const mesh::ReconstructPoint& interpolation = metrics->copy_pt(f); // <------ double check!!!
              neighbor_flx_pt_vars[flx_pt].setZero();
              boost_foreach( Uint n, interpolation.used_points() )
              {
                const Real C = interpolation.coeff(n);
                for (Uint v=0; v<nb_vars; ++v)
                {
                  neighbor_flx_pt_vars[flx_pt][v] += C * m_field->array()[nodes[n]][v];
                }
              }



              const std::vector<mesh::ReconstructPoint>& gradient = metrics->zero_grad(); // <------ double check!!!
              neighbor_flx_pt_grad[flx_pt].setZero();
              for (Uint d=0; d<NDIM; ++d)
              {
                boost_foreach( Uint n, gradient[d].used_points() )
                {
                  const Real C = gradient[d].coeff(n);
                  for (Uint v=0; v<nb_vars; ++v)
                  {
                    neighbor_flx_pt_grad[flx_pt](d,v) += C * m_field->array()[nodes[n]][v];
                  }
                }
              }

              // ---> decide to just use inside gradient instead!
              neighbor_flx_pt_grad[flx_pt] = flx_pt_grad[flx_pt];
//              neighbor_flx_pt_grad[flx_pt] = metrics->flx_pt_Jinv(flx_pt)*neighbor_flx_pt_grad[flx_pt];
//              CFinfo << neighbor_flx_pt_grad[flx_pt].transpose() << CFendl;

//              m_term->get_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
//                                     /*in*/  m_metrics->flx_pt_coords(flx_pt),
//                                     /*in*/  m_metrics->copy_pt(f),
//                                     /*in*/  m_metrics->zero_grad(),        // <----double check!!!
//                                     /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
//                                     /*out*/ m_neighbor_flx_pt_vars[flx_pt],
//                                     /*out*/ m_neighbor_flx_pt_gradvars[flx_pt],
//                                     /*out*/ m_neighbor_flx_pt_gradvars_grad[flx_pt] );

//              m_neighbor_flx_pt_gradvars_grad[flx_pt] = m_flx_pt_gradvars_grad[flx_pt];
            }
            else // Gradient is computed in neighbor cell
            {
              mesh::Connectivity::ConstRow nodes = neighbor_space[face_nb]->connectivity()[neighbor_elem_idx[face_nb]];

              const mesh::ReconstructPoint& interpolation = metrics->interpolation_from_sol_pts_to_flx_pt(neighbor_flx_pt);
              neighbor_flx_pt_vars[flx_pt].setZero();
              boost_foreach( Uint n, interpolation.used_points() )
              {
                const Real C = interpolation.coeff(n);
                for (Uint v=0; v<nb_vars; ++v)
                {
                  neighbor_flx_pt_vars[flx_pt][v] += C * m_field->array()[nodes[n]][v];
                }
              }



              const std::vector<mesh::ReconstructPoint>& gradient = metrics->gradient_from_sol_pts_to_flx_pt(neighbor_flx_pt);
              neighbor_flx_pt_grad[flx_pt].setZero();
              for (Uint d=0; d<NDIM; ++d)
              {
                boost_foreach( Uint n, gradient[d].used_points() )
                {
                  const Real C = gradient[d].coeff(n);
                  for (Uint v=0; v<nb_vars; ++v)
                  {
                    neighbor_flx_pt_grad[flx_pt](d,v) += C * m_field->array()[nodes[n]][v];
                  }
                }
              }
//              neighbor_flx_pt_grad[flx_pt] = metrics->flx_pt_Jinv(flx_pt)*neighbor_flx_pt_grad[flx_pt];
//              CFinfo << neighbor_flx_pt_grad[flx_pt].transpose() << CFendl;


//              m_term->get_variables( /*in*/  *m_neighbor_space[face_nb], m_neighbor_elem_idx[face_nb],
//                                     /*in*/  m_metrics->flx_pt_coords(flx_pt),
//                                     /*in*/  m_metrics->interpolation_from_sol_pts_to_flx_pt(neighbor_flx_pt),
//                                     /*in*/  m_metrics->gradient_from_sol_pts_to_flx_pt(neighbor_flx_pt),
//                                     /*in*/  m_metrics->flx_pt_Jinv(flx_pt),
//                                     /*out*/ m_neighbor_flx_pt_vars[flx_pt],
//                                     /*out*/ m_neighbor_flx_pt_gradvars[flx_pt],
//                                     /*out*/ m_neighbor_flx_pt_gradvars_grad[flx_pt] );
            }

            avg_flx_pt_vars[flx_pt] += neighbor_flx_pt_vars[flx_pt];
            avg_flx_pt_vars[flx_pt] *= 0.5;

          } // end for (f) face_pts


          // Compute Left lifting operator ( Lambda_L )
          boost_foreach(RealVector& j, jump)
          {
            j.setZero();
          }
          boost_foreach(RealVector& j, neighbor_jump)
          {
            j.setZero();
          }


          for (Uint f=0; f<nb_face_pts; ++f)
          {
            Uint flx_pt = face_pts[face_nb][f];
            jump[flx_pt] = neighbor_flx_pt_vars[flx_pt] - flx_pt_vars[flx_pt];
          }

          for (Uint f=0; f<nb_face_pts; ++f)
          {
            Uint flx_pt = face_pts[face_nb][f];
            LambdaL[flx_pt].setZero();
            for (Uint d=0; d<NDIM; ++d)
            {
              boost_foreach( Uint n, metrics->derivation_from_flx_pts_to_flx_pt(flx_pt,d).used_points())
              {
                const Real C = metrics->derivation_from_flx_pts_to_flx_pt(flx_pt,d).coeff(n);
                for (Uint v=0; v<nb_vars; ++v)
                  LambdaL[flx_pt](d,v) += jump[n][v] * C;
              }
            }
//            LambdaL[flx_pt] = metrics->flx_pt_Jinv(flx_pt) * LambdaL[flx_pt];
          }

          // Compute Right lifting operator ( Lambda_R )
          if ( connected->is_bdry_face()[face_nb] )
          {
            for (Uint f=0; f<nb_face_pts; ++f)
            {
              Uint flx_pt = face_pts[face_nb][f];
              LambdaR[flx_pt] = LambdaL[flx_pt];
            }
          }
          else
          {
            for (Uint f=0; f<nb_face_pts; ++f)
            {
              Uint flx_pt = face_pts[face_nb][f];
              Uint neighbor_flx_pt = neighbor_face_pts[face_nb][f];
              neighbor_jump[neighbor_flx_pt] = - jump[flx_pt];
            }
            for (Uint f=0; f<nb_face_pts; ++f)
            {
              Uint flx_pt = face_pts[face_nb][f];
              Uint neighbor_flx_pt = neighbor_face_pts[face_nb][f];
              LambdaR[flx_pt].setZero();
              for (Uint d=0; d<NDIM; ++d)
              {
                boost_foreach( Uint n, metrics->derivation_from_flx_pts_to_flx_pt(neighbor_flx_pt,d).used_points() )
                {
                  const Real C = metrics->derivation_from_flx_pts_to_flx_pt(neighbor_flx_pt,d).coeff(n);
                  for (Uint v=0; v<nb_vars; ++v)
                    LambdaR[flx_pt](d,v) += neighbor_jump[n][v] * C;
                }
              }
//              LambdaR[flx_pt] = metrics->flx_pt_Jinv(flx_pt) * LambdaR[flx_pt];
            }
          }
          // Now LambdaL and LambdaR are computed

          // Average now
          for (Uint f=0; f<nb_face_pts; ++f) // loop over face flux points
          {
            Uint flx_pt = face_pts[face_nb][f];
            avg_flx_pt_grad[flx_pt] += neighbor_flx_pt_grad[flx_pt];
            avg_flx_pt_grad[flx_pt] += m_alpha*LambdaL[flx_pt];
            avg_flx_pt_grad[flx_pt] += m_alpha*LambdaR[flx_pt];
            avg_flx_pt_grad[flx_pt] *= 0.5;
          }
        } // end for (face_nb)

        std::vector<Eigen::Matrix<Real,NDIM,Eigen::Dynamic> > grad_pt_grad(nb_grad_pts, Eigen::Matrix<Real,NDIM,Eigen::Dynamic>(NDIM,nb_vars));
        for (Uint grad_pt=0; grad_pt<nb_grad_pts; ++grad_pt)
        {
          grad_pt_grad[grad_pt].setZero();
          for (Uint d=0; d<NDIM; ++d)
          {
            boost_foreach( Uint flx_pt, line_interpolation_from_flx_pts_to_grad_pt[grad_pt][d].used_points() )
            {
              const Real C = line_interpolation_from_flx_pts_to_grad_pt[grad_pt][d].coeff(flx_pt);
              for (Uint v=0; v<nb_vars; ++v)
                grad_pt_grad[grad_pt](d,v) += avg_flx_pt_grad[flx_pt](d,v) * C;
            }
          }
          grad_pt_grad[grad_pt] = cells->element_type().jacobian(grad_sf->local_coordinates().row(grad_pt), metrics->cell_coords()).inverse() * grad_pt_grad[grad_pt];
        }

        mesh::Connectivity::ConstRow grad_nodes = grad_space->connectivity()[e];
        for (Uint grad_pt=0; grad_pt<nb_grad_pts; ++grad_pt)
        {
          Uint p = grad_nodes[grad_pt];
          for (Uint d=0; d<NDIM; ++d)
          {
            for (Uint v=0; v<nb_vars; ++v)
            {
              m_field_gradient->array()[p][v+d*nb_vars] = grad_pt_grad[grad_pt](d,v);
            }
          }
        }


      } // end for (e)
    } // end if cells
  } // end for (space)
}

////////////////////////////////////////////////////////////////////////////////

} // tools
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_tools_ComputeFieldGradientBR2_hpp
