// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_test_Term_hpp
#define cf3_sdm_test_Term_hpp

////////////////////////////////////////////////////////////////////////////////

#include <boost/bind.hpp>
#include "cf3/common/Component.hpp"
#include "cf3/solver/Term.hpp"

namespace cf3 {
namespace sdm {
namespace test {

////////////////////////////////////////////////////////////////////////////////

template < Uint NB_DIM >
class Term : public solver::TermBase<NB_DIM,NB_DIM,NB_DIM,NB_DIM>
{
public: 

  /// @brief Constructor
  Term( const std::string& name ) : solver::TermBase<NB_DIM,NB_DIM,NB_DIM,NB_DIM>(name) {}
  
  /// @brief Destructor
  virtual ~Term() {}
  
  static std::string type_name() { return "Term"; };
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_DIFFUSION  = true };
  enum { ENABLE_SOURCE     = true };
  

  struct DATA
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM coords;
  };


  /// @brief Compute variables and gradients in a given element point
  ///
  /// The interpolation and gradient reconstructions, as well as
  void get_variables( const mesh::Space& space,
                      const Uint elem_idx,
                      const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& coords,
                      const mesh::ReconstructPoint& interpolation,
                      const std::vector<mesh::ReconstructPoint>& gradient,
                      const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNDIM& jacobian,
                      const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNDIM& jacobian_inverse,
                      const Real& jacobian_determinant,
                      typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NVAR& vars,
                      typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NGRAD& gradvars,
                      typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNGRAD& gradvars_grad )
  {
    // Typically get the solution here
    RealMatrix cell_coords = space.get_coordinates(elem_idx);
    vars.setZero();
    boost_foreach( Uint pt, interpolation.used_points() )
    {
      for (Uint v=0; v<NB_DIM; ++v)
      {
        vars[v] += interpolation.coeff(pt) * cell_coords(pt,v);
      }
    }

    // because NGRAD = NVAR in this case, and to make a point
    gradvars = vars;
    gradvars_grad.setZero();
    for (Uint d=0; d<NB_DIM; ++d)
    {
      boost_foreach( Uint pt, gradient[d].used_points() )
      {
        for (Uint v=0; v<NB_DIM; ++v)
        {
          gradvars_grad(d,v) += gradient[d].coeff(pt) * cell_coords(pt,v);
        }
      }
    }
    gradvars_grad = jacobian_inverse * gradvars_grad;
  }

  /// @brief Compute bdry variables and gradients in a given element point
  void get_bdry_variables( const mesh::Space& space,
                           const Uint elem_idx,
                           const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& coords,
                           const mesh::ReconstructPoint& interpolation,
                           const std::vector<mesh::ReconstructPoint>& gradient,
                           const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNDIM& jacobian,
                           const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNDIM& jacobian_inverse,
                           const Real& jacobian_determinant,
                           typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NVAR& vars,
                           typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NGRAD& gradvars,
                           typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNGRAD& gradvars_grad )
  {
    // Typically get the solution here
    RealMatrix cell_coords = space.get_coordinates(elem_idx);
    vars.setZero();
    boost_foreach( Uint pt, interpolation.used_points() )
    {
      for (Uint v=0; v<NB_DIM; ++v)
      {
        vars[v] += interpolation.coeff(pt) * cell_coords(pt,v);
      }
    }
    vars = coords;

    // because NGRAD = NVAR in this case, and to make a point
    gradvars = vars;
    gradvars_grad.setZero();
    for (Uint d=0; d<NB_DIM; ++d)
    {
      gradvars_grad(d,d)=1.;
    }
  }

  
  /// @brief Set constants in the data
  void set_phys_data_constants( DATA& phys_data ) { }

  /// @brief Compute the data from computed variables and gradients
  void compute_phys_data( const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& coords,
                          const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NVAR& vars,
                          const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NGRAD& gradvars,
                          const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data )
  {
    static const Real eps = 4000*math::Consts::eps();
    for (Uint d=0; d<NB_DIM; ++d)
    {
      if ( std::abs(vars[d]-coords[d]) > eps )
      {
        throw common::BadValue( FromHere(), "Interpolations or face-connectivities failed (difference > "+common::to_str(eps)+")");
      }
      for (Uint v=0; v<NB_DIM; ++v)
      {
        if (v == d)
        {
          if ( std::abs(gradvars_grad(d,v)-1.) > eps )
          {
            std::cout << gradvars_grad << std::endl;
            throw common::BadValue( FromHere(), "Gradient computation failed (difference > "+common::to_str(eps)+")");
          }
        }
        else
        {
          if ( std::abs(gradvars_grad(d,v)-0.) > eps )
          {
            std::cout << gradvars_grad << std::endl;
            throw common::BadValue( FromHere(), "Gradient computation failed (difference > "+common::to_str(eps)+")");
          }
        }
      }
    }
    phys_data.coords = vars;
  }
  
public: // flux computations

  static void compute_convective_flux( const DATA& p, const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& normal,
                                       typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NEQS& flux, Real& wave_speed )
  {
    flux.setZero();
    wave_speed=0.;
  }

  static  void compute_riemann_flux( const DATA& left, const DATA& right, const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& normal,
                                     typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NEQS& flux, Real& wave_speed )
  {
    flux.setZero();
    wave_speed=0.;
  }
  
  static void compute_diffusive_flux( const DATA& p, const typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::ColVector_NDIM& normal,
                                      typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NEQS& flux, Real& wave_speed )
  {
    flux.setZero();
    wave_speed=0.;
  }
  
  static void compute_source( const DATA& p, typename physics::MatrixTypes<NB_DIM,NB_DIM,NB_DIM,NB_DIM>::RowVector_NEQS& source )
  {
    source = p.coords;
  }
};

////////////////////////////////////////////////////////////////////////////////

} // test
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_test_Term_hpp
