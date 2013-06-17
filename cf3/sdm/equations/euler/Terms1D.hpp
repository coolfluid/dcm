// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_euler_Terms1D_hpp
#define cf3_sdm_equations_euler_Terms1D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/math/Checks.hpp"
#include "cf3/solver/Term.hpp"
#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/sdm/equations/euler/LibEuler.hpp"
#include "cf3/physics/euler/euler1d/Functions.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API Terms1D : public solver::TermBase< 1 /*dim*/, 3 /*eqs*/, 3 /*vars*/, 0/*grads*/ >
{
public: 

  /// @brief Constructor
  Terms1D( const std::string& name );
  
  /// @brief Destructor
  virtual ~Terms1D() {}
  
  static std::string type_name() { return "Terms1D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };

  typedef physics::euler::euler1d::Data DATA;

public: // Variable and PhysData computation
    
  /// @brief Compute variables and gradients in a given element point
  ///
  /// The interpolation and gradient reconstructions, as well as
  void get_variables( const mesh::Space& space,
                      const Uint elem_idx,
                      const ColVector_NDIM& coords,
                      const mesh::ReconstructPoint& interpolation,
                      const std::vector<mesh::ReconstructPoint>& gradient,
                      const Matrix_NDIMxNDIM& jacobian,
                      const Matrix_NDIMxNDIM& jacobian_inverse,
                      const Real& jacobian_determinant,
                      RowVector_NVAR& vars,
                      RowVector_NGRAD& gradvars,
                      Matrix_NDIMxNGRAD& gradvars_grad );

  void get_bdry_variables( const mesh::Space& space,
                           const Uint elem_idx,
                           const ColVector_NDIM& coords,
                           const mesh::ReconstructPoint& interpolation,
                           const std::vector<mesh::ReconstructPoint>& gradient,
                           const Matrix_NDIMxNDIM& jacobian,
                           const Matrix_NDIMxNDIM& jacobian_inverse,
                           const Real& jacobian_determinant,
                           RowVector_NVAR& vars,
                           RowVector_NGRAD& gradvars,
                           Matrix_NDIMxNGRAD& gradvars_grad );
  
  /// @brief Set constants in the data
  void set_phys_data_constants( DATA& phys_data );

  /// @brief Compute the data from computed variables and gradients
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Flux computations

  static void compute_convective_flux( const DATA& p, const ColVector_NDIM& normal,
                                       RowVector_NEQS& flux, Real& wave_speed )
  {   
    physics::euler::euler1d::compute_convective_flux(p,normal,flux,wave_speed);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      if (math::Checks::is_nan(flux[eq]))
      {
        std::cout << "flux = " << flux << std::endl;
        throw common::FailedToConverge(FromHere(), "flux is nan");
      }
    }
  }

  void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                             RowVector_NEQS& flux, Real& wave_speed )
  {
    m_riemann_solver->compute_riemann_flux(left,right,normal,flux,wave_speed);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      if (math::Checks::is_nan(flux[eq]))
      {
        std::cout << "QL = " << left.cons << std::endl;
        std::cout << "QR = " << right.cons << std::endl;

        std::cout << "flux = " << flux << std::endl;
        throw common::FailedToConverge(FromHere(), "flux is nan");
      }
    }
  }

private:

  void config_riemann_solver();

private: // configuration

  Handle< solver::RiemannSolver<Terms1D> > m_riemann_solver;
  std::string m_riemann_solver_type;
  Real m_gamma;
  Real m_R;
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_euler_Terms1D_hpp
