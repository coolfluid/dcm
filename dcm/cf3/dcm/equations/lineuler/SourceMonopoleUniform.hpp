// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_SourceMonopoleUniform2D_hpp
#define cf3_dcm_equations_lineuler_SourceMonopoleUniform2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/Time.hpp"
#include "cf3/solver/Term.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {


////////////////////////////////////////////////////////////////////////////////

/// @brief Monopole source term for the 2D Linearized Euler Equations
///        with uniform mean flow
/// @author Willem Deconinck
///
/// @f[ f(x,y,t) = \sin(\omega\ t) \  \varepsilon\ \exp(-\alpha\ ( (x-x0)^2 + y-y0)^2 ) ) @f]
/// @f[ S = \[ \begin{array}[c]
///     f(x,y,t) \\
///     0        \\
///     0        \\
///     c_0^2 \ f(x,y,t)
/// \end{array} \] @f]
///
/// Heat capacity ratio gamma and mean flow is necessary to compute the sound speed c0
class dcm_equations_lineuler_API SourceMonopoleUniform2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 0 /*vars*/, 0/*grads*/ >
{
public: 

  /// @brief Constructor
  SourceMonopoleUniform2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~SourceMonopoleUniform2D() {}
  
  static std::string type_name() { return "SourceMonopoleUniform2D"; }

public: // types

  enum { ENABLE_SOURCE = true };

  struct DATA
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ColVector_NDIM coords;
    Real c02;
  };

public: // Variable and PhysData computation
  
  /// @brief Set constants in the data
  void set_phys_data_constants( DATA& phys_data );
  
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Source computation

  void compute_source( const DATA& p, RowVector_NEQS& source );
  
private: // functions

  /// spatial part of the source term
  Real f(const ColVector_NDIM& coord);

private:
  
  Real m_source;                      ///< dummy variable

  Real m_alpha;                       ///< Width of monopole
  Real m_eps;                         ///< Amplitude of monopole
  Real m_freq;                        ///< Frequency of monopole
  std::vector<Real> m_source_loc;     ///< Location of monopole

  Real m_gamma;                       ///< Heat capacity ratio
  Real m_rho0;
  Real m_p0;
};

////////////////////////////////////////////////////////////////////////////////

/// @brief Monopole source term for the 3D Linearized Euler Equations
///        with uniform mean flow
/// @author Willem Deconinck
///
/// @f[ f(x,y,z,t) = \sin(\omega\ t) \  \varepsilon\ \exp(-\alpha\ ( (x-x0)^2 + y-y0 + z-z0)^2 ) ) @f]
/// @f[ S = \[ \begin{array}[c]
///     f(x,y,t) \\
///     0        \\
///     0        \\
///     0        \\
///     c_0^2 \ f(x,y,z,t)
/// \end{array} \] @f]
///
/// Heat capacity ratio gamma and mean flow is necessary to compute the sound speed c0
class dcm_equations_lineuler_API SourceMonopoleUniform3D : public solver::TermBase< 3 /*dim*/, 5 /*eqs*/, 0 /*vars*/, 0/*grads*/ >
{
public: 

  /// @brief Constructor
  SourceMonopoleUniform3D( const std::string& name );
  
  /// @brief Destructor
  virtual ~SourceMonopoleUniform3D() {}
  
  static std::string type_name() { return "SourceMonopoleUniform3D"; }
  
public: // types

  enum { ENABLE_SOURCE = true };

  struct DATA
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ColVector_NDIM coords;
    Real c02;
  };
  
public: // Variable and PhysData computation
  
  /// @brief Set constants in the data
  void set_phys_data_constants( DATA& phys_data );
  
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Source computation

  void compute_source( const DATA& p, RowVector_NEQS& source );
  
private: // functions

  /// spatial part of the source term
  Real f(const ColVector_NDIM& coord);

private:
  
  Real m_source;                      ///< dummy variable

  Real m_alpha;                       ///< Width of monopole
  Real m_eps;                         ///< Amplitude of monopole
  Real m_freq;                        ///< Frequency of monopole
  std::vector<Real> m_source_loc;     ///< Location of monopole

  Real m_gamma;                       ///< Heat capacity ratio
  Real m_rho0;
  Real m_p0;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_SourceMonopoleUniform3D_hpp
