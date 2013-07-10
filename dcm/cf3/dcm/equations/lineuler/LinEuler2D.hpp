// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_LinEuler2D_hpp
#define cf3_dcm_equations_lineuler_LinEuler2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/Physics.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"
#include "cf3/dcm/core/Reconstructions.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

/// @brief Linearized Euler physics, for Spectral Difference method
///
/// @author Willem Deconinck
class dcm_equations_lineuler_API LinEuler2D : public dcm::Physics {

public:

  enum {NDIM = 2};
  enum {NEQS = 4};

  typedef Eigen::Matrix<Real,NDIM,1> RealVectorNDIM;
  typedef Eigen::Matrix<Real,NEQS,1> RealVectorNEQS;
  typedef Eigen::Matrix<Real,NDIM,NDIM,Eigen::RowMajor> RealMatrixNDIMxNDIM;
  typedef Eigen::Matrix<Real,NEQS,NDIM,Eigen::RowMajor> RealMatrixNEQSxNDIM;
  typedef Eigen::Matrix<Real,NDIM,NEQS,Eigen::RowMajor> RealMatrixNDIMxNEQS;
  typedef Eigen::Matrix<Real,NEQS,NEQS,Eigen::RowMajor> RealMatrixNEQSxNEQS;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  struct PhysData
  {
    // To be set manually:
    LinEuler2D::RealVectorNDIM coord;
    LinEuler2D::RealVectorNEQS solution;

    Real rho0;
    LinEuler2D::RealVectorNDIM U0;
    Real p0;

    Real gamma;

//    RealMatrix2 grad_U0;   // [ du0/dx , dv0/dx ; du0/dy , dv0/dy ]
//    RealVector2 grad_p0;   // [ dp0/dx ; dp0/dy ]
    void compute()
    {
      rho = solution[0];
      rho0_U[XX] = solution[1];
      rho0_U[YY] = solution[2];
      p = solution[3];
      U[XX] = rho0_U[XX]/rho0;
      U[YY] = rho0_U[YY]/rho0;

      c0 = std::sqrt(gamma*p0/rho0);
    }

    // To be computed
    Real c0;
    Real rho;
    LinEuler2D::RealVectorNDIM rho0_U;
    Real p;
    LinEuler2D::RealVectorNDIM U;

//    RealVector2 grad_rho;      // [ drho/dx ; drho/dy ]
//    RealMatrix2 grad_rho0_U;   // [ drho0_u/dx , drho0_v/dx ; drho0_u/dy , drho0_v/dy ]
//    RealVector2 grad_p;        // [ dp/dx ; dp/dy ]
  };

public: // functions

  /// Contructor
  /// @param name of the component
  LinEuler2D ( const std::string& name );

  /// Virtual destructor
  virtual ~LinEuler2D();

  /// Get the class name
  static std::string type_name () { return "LinEuler2D"; }
  
  virtual void create_terms( const Handle<dcm::DomainDiscretization>& domain_discretization );
  
  virtual void create_fields( const Handle<mesh::Dictionary>& dict );

  void compute_sol_pt_phys_data( const mesh::SpaceElem& cell, const Uint sol_pt,
                                 PhysData& phys_data );

  void compute_flx_pt_phys_data( const mesh::SpaceElem& cell, const Uint flx_pt,
                                 PhysData& phys_data );

  static void compute_interior_fluxes(const PhysData& phys_data, const RealVectorNDIM& normal,
                                      RealVectorNEQS& flux,
                                      Real& wave_speed);

  static void compute_face_fluxes(const PhysData& left, const PhysData& right, const RealVectorNDIM &normal,
                                   RealVectorNEQS& flux,
                                   Real& wave_speed);

  void init_fields();

  static void compute_convective_flux( const PhysData& phys_data, const RealVectorNDIM& normal,
                                       RealVectorNEQS& flux);

  static void compute_convective_wave_speed( const PhysData& phys_data, const RealVectorNDIM& normal,
                                             Real& wave_speed);

  static void compute_convective_eigenvalues( const PhysData& phys_data, const RealVectorNDIM& normal,
                                              RealVectorNEQS& eigenvalues);

  static void compute_convective_eigensystem( const PhysData& p, const RealVectorNDIM& normal,
                                              RealMatrixNEQSxNEQS& right_eigenvectors,
                                              RealVectorNEQS& eigenvalues,
                                              RealMatrixNEQSxNEQS& left_eigenvectors);

  /// @brief Optimized version to compute the absolute flux jacobian
  ///
  /// This is possible because only the mean flow variables are present in this expression,
  /// and thus no averaged state needs to be constructed to evaluate this.
  /// Equivalent to  right_eigenvectors*abs(eigenvalues)*left_eigenvectors
  static void compute_absolute_flux_jacobian( const PhysData& p, const RealVectorNDIM& normal,
                                              RealMatrixNEQSxNEQS& absolute_flux_jacobian);

private: // functions

  void configure_terms();




private: // data

  Handle<dcm::Term>                         m_convection;
  Handle<dcm::Term>                         m_mean_flow_source;
  Handle<solver::Time>                      m_time;

  Handle<mesh::Dictionary>                  m_solution_space;

  Handle<mesh::Field>                       m_coords;
  Handle<mesh::Field>                       m_solution;
  Handle<mesh::Field>                       m_rho0;
  Handle<mesh::Field>                       m_U0;
  Handle<mesh::Field>                       m_p0;
  Handle<mesh::Field>                       m_grad_U0;
  Handle<mesh::Field>                       m_grad_p0;

  Real m_gamma;

  std::map<mesh::Space*,ReconstructToFluxPoints> m_solution_to_flux;
  std::map<mesh::Space*,ReconstructToFluxPoints> m_geometry_to_flux;

  // Variables to be used in static functions
  static RealMatrixNEQSxNEQS A;
  static RealMatrixNEQSxNEQS L;
  static RealMatrixNEQSxNEQS R;
  static RealVectorNEQS D;
  static RealVectorNEQS flux_left;
  static RealVectorNEQS flux_right;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_LinEuler2D_hpp
