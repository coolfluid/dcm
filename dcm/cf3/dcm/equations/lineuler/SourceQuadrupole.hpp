// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_SourceQuadrupole2D_hpp
#define cf3_dcm_equations_lineuler_SourceQuadrupole2D_hpp

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

/// @brief Quadrupole source term for the 2D Linearized Euler Equations
///        with uniform mean flow
/// @author Willem Deconinck
///
///  @f[ S = \[ \begin{array}[c]
///     0 \\
///      \epsilon*\sin((\pi/4*(x-x_0)/\alpha) * \exp(-\log(2)/\alpha*(y-y_0)^2) * \sin(\omega*t) \\
///     -\epsilon*\sin((\pi/4*(y-y_0)/\alpha) * \exp(-\log(2)/\alpha*(x-x_0)^2) * \sin(\omega*t) \\
///     0
/// \end{array} \] @f]
class dcm_equations_lineuler_API SourceQuadrupole2D : public solver::TermBase< 2 /*dim*/, 4 /*eqs*/, 0 /*vars*/, 0/*grads*/ >
{
public: 

  /// @brief Constructor
  SourceQuadrupole2D( const std::string& name );
  
  /// @brief Destructor
  virtual ~SourceQuadrupole2D() {}
  
  static std::string type_name() { return "SourceQuadrupole2D"; }

public: // types

  enum { ENABLE_SOURCE = true };

  struct DATA
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ColVector_NDIM coords;
  };

public: // Variable and PhysData computation
  
  void compute_phys_data( const ColVector_NDIM& coords,
                          const RowVector_NVAR& vars,
                          const RowVector_NGRAD& gradvars,
                          const Matrix_NDIMxNGRAD& gradvars_grad,
                          DATA& phys_data );
  
public: // Source computation

  void compute_source( const DATA& p, RowVector_NEQS& source );
  
private: // functions

  /// spatial part of the source term
  Real fx(const ColVector_NDIM& coord);

  /// spatial part of the source term
  Real fy(const ColVector_NDIM& coord);

private:

  Real m_theta;                       ///< Angle of directivity
  Real m_alpha;                       ///< Width of quadrupole
  Real m_eps;                         ///< Amplitude of quadrupole
  Real m_freq;                        ///< Frequency of quadrupole
  std::vector<Real> m_source_loc;     ///< Location of quadrupole
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_SourceQuadrupole_hpp
