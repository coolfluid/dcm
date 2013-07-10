// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_navierstokes_RightHandSide2D_hpp
#define cf3_dcm_equations_navierstokes_RightHandSide2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/Term.hpp"
#include "cf3/solver/RiemannSolver.hpp"
#include "cf3/physics/navierstokes/navierstokes2d/Functions.hpp"
#include "cf3/dcm/equations/navierstokes/LibNavierStokes.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_navierstokes_API RightHandSide2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 4 /*vars*/, 3/*grads*/ >
{
public: 

  /// @brief Constructor
  RightHandSide2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~RightHandSide2D() {}
  
  static std::string type_name() { return "NSRightHandSide2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true };
  enum { ENABLE_DIFFUSION  = true };

  typedef physics::navierstokes::navierstokes2d::Data DATA;

public: // Variable and PhysData computation

  // About 50% of execution is spent in this class for P1 with LUSGS

  /// @brief Compute variables and gradients in a given element point
  // 32% of execution!!! WHY
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
  //  12% of execution!!!
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Flux computations

  // 0.2% execution
  static void compute_convective_flux( const DATA& p, const ColVector_NDIM& normal,
                                       RowVector_NEQS& flux, Real& wave_speed )
  {
    physics::euler::euler2d::compute_convective_flux(p,normal,flux,wave_speed);
  }

  // 5.6% execution
  void compute_riemann_flux( const DATA& left, const DATA& right, const ColVector_NDIM& normal,
                             RowVector_NEQS& flux, Real& wave_speed )
  {
    m_riemann_solver->compute_riemann_flux(left,right,normal,flux,wave_speed);
  }

  // 1.5% execution
  static void compute_diffusive_flux(const DATA &p, const ColVector_NDIM &normal,
                                     RowVector_NEQS &flux, Real &wave_speed);


private:

  void config_riemann_solver();

private: // configuration

  Handle< solver::RiemannSolver<cf3::physics::euler::euler2d::Data,NDIM,NEQS> > m_riemann_solver;
  std::string m_riemann_solver_type;
  Real m_gamma;
  Real m_R;
  Real m_kappa;
  Real m_mu;

  Real _rho;
  Real _u;
  Real _v;
  Real _E;
  Real _U2;
  Real _p;
  Real _T;
  ColVector_NDIM _grad_rho;
  ColVector_NDIM _grad_rhou;
  ColVector_NDIM _grad_rhov;
  ColVector_NDIM _grad_rhoE;
  ColVector_NDIM _grad_u;
  ColVector_NDIM _grad_v;
  ColVector_NDIM _grad_T;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_navierstokes_API Convection2D : public RightHandSide2D
{
public: 

  /// @brief Constructor
  Convection2D( const std::string& name ) : RightHandSide2D(name) {}
  
  /// @brief Destructor
  virtual ~Convection2D() {}
  
  static std::string type_name() { return "NSConvection2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = true  };
  enum { ENABLE_DIFFUSION  = false };

};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_navierstokes_API Diffusion2D : public RightHandSide2D
{
public: 

  /// @brief Constructor
  Diffusion2D( const std::string& name ) : RightHandSide2D(name) {}
  
  /// @brief Destructor
  virtual ~Diffusion2D() {}
  
  static std::string type_name() { return "NSDiffusion2D"; }
  
public: // types

  enum { ENABLE_CONVECTION = false };
  enum { ENABLE_DIFFUSION  = true  };

};

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_navierstokes_RightHandSide2D_hpp
