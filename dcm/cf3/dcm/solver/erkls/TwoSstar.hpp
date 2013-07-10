// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_solver_erkls_TwoSstar_hpp
#define cf3_dcm_solver_erkls_TwoSstar_hpp

#include "cf3/solver/PDESolver.hpp"
#include "cf3/dcm/solver/erkls/LibERKLS.hpp"
#include "cf3/dcm/solver/erkls/TwoSstarCoeffs.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace solver {
namespace erkls {

/////////////////////////////////////////////////////////////////////////////////////

/// Runge-Kutta low storage integration method using only 2 registers
/// @ref David I. Ketcheson: Runge-Kutta methods with minimum storage implementations
///      Journal of Computational Physics 229 (2010) 1763–1773
///      doi:10.1016/j.jcp.2009.11.006
/// The order is not necessarily the same as the number of stages "m"
/// The order depends on the coefficients alpha and beta
/// Algorithm 2S* with m = number of stages (not necessarily same as order)
/// @code
/// S1 := U(t=n)   S2 := U(t=n)
/// for i = 2:m+1 do
///    S1 := (1-alpha(i,1))*S1 + alpha(i,1)*S2 + beta(i,i-1)*dt*F(S1)
/// end
/// U(t=n+1) = S1
/// @endcode
/// @author Willem Deconinck
class dcm_solver_erkls_API TwoSstarImplementation : public cf3::solver::PDESolver {

public: // functions

  /// Contructor
  /// @param name of the component
  TwoSstarImplementation ( const std::string& name );

  /// Virtual destructor
  virtual ~TwoSstarImplementation() {}

  /// Get the class name
  static std::string type_name () { return "TwoSstarImplementation"; }

  /// execute the action
  virtual void step ();

private: // functions

  void create_fields();

protected: // data

  Handle<TwoSstarCoeffs> m_coeffs;

private: // data

  /// Second register necessary for low-storage runge kutta algorithm  2S*
  Handle<mesh::Field> m_backup; // ( = S2 in the above @code )
  Handle<mesh::Field> m_dt;
};

/////////////////////////////////////////////////////////////////////////////////////

class dcm_solver_erkls_API TwoSstar : public TwoSstarImplementation {

public: // functions

  /// Contructor
  /// @param name of the component
  TwoSstar ( const std::string& name );

  /// Virtual destructor
  virtual ~TwoSstar() {}

  /// Get the class name
  static std::string type_name () { return "TwoSstar"; }

private: // functions

  void config_order();

};

/////////////////////////////////////////////////////////////////////////////////////

template <typename COEFFS>
class TwoSstarT : public TwoSstarImplementation {
public:

  /// @brief Type name
  static std::string type_name () { return COEFFS::name(); }

  /// @brief Constructor
  TwoSstarT(const std::string& name) : TwoSstarImplementation(name)
  {
    m_coeffs = create_component<TwoSstarCoeffs>("coeffs");
    m_coeffs->set( COEFFS() );
  }
};

/////////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // dcm
} // cf3

#endif // cf3_dcm_solver_erkls_TwoSstar_hpp
