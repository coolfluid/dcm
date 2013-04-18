// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_solver_erkls_ThreeSstar_hpp
#define cf3_sdm_solver_erkls_ThreeSstar_hpp

#include "cf3/common/OptionList.hpp"
#include "cf3/solver/PDESolver.hpp"
#include "cf3/solver/TimeStepComputer.hpp"
#include "cf3/sdm/solver/erkls/LibERKLS.hpp"
#include "cf3/sdm/solver/erkls/ThreeSstarCoeffs.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace solver {
namespace erkls {

/////////////////////////////////////////////////////////////////////////////////////

/// Runge-Kutta low storage integration method using only 3 registers
/// @ref David I. Ketcheson: Runge-Kutta methods with minimum storage
///      implementations. Journal of Computational Physics 229 (2010) 1763–1773
///      doi:10.1016/j.jcp.2009.11.006
/// The order is not necessarily the same as the number of stages "m"
/// The order depends on the coefficients alpha and beta
/// Algorithm 3S* with m = number of stages (not necessarily same as order)
/// @code
/// S1 := U(t=n)   S2 := 0   S3 := U(t=n)
/// for i = 2:m+1 do
///     S2 := S2 + delta(i-1)*S1
///     S1 := gamma(i,1)*S1 + gamma(i,2)*S2 + gamma(i,3)*S3 + beta(i,i-1)*dt*F(S1)
/// end
/// U(t=n+1) = S1
/// @endcode
/// @author Willem Deconinck
/// @author Matteo Parsani
class sdm_solver_erkls_API ThreeSstarImplementation : public cf3::solver::PDESolver {

public: // functions

  /// Contructor
  /// @param name of the component
  ThreeSstarImplementation ( const std::string& name );

  /// Virtual destructor
  virtual ~ThreeSstarImplementation() {}

  /// Get the class name
  static std::string type_name () { return "ThreeSstarImplementation"; }

  /// execute the action
  virtual void step ();

private: // functions

  void create_fields();

protected: // data

  Handle<ThreeSstarCoeffs> m_coeffs;

private: // data

  /// Second register necessary for low-storage runge kutta algorithm 3S*
  Handle<mesh::Field> m_S2;
  /// Third register necessary for low-storage runge kutta algorithm  3S*
  Handle<mesh::Field> m_backup; // ( = S3 in algorithm )
  Handle<mesh::Field> m_dt;
};

/////////////////////////////////////////////////////////////////////////////////////

class sdm_solver_erkls_API ThreeSstar : public ThreeSstarImplementation {

public: // functions

  /// Contructor
  /// @param name of the component
  ThreeSstar ( const std::string& name );

  /// Virtual destructor
  virtual ~ThreeSstar() {}

  /// Get the class name
  static std::string type_name () { return "ThreeSstar"; }

};

/////////////////////////////////////////////////////////////////////////////////////

template <typename COEFFS>
class ThreeSstarT : public ThreeSstarImplementation {
public:

  /// @brief Type name
  static std::string type_name () { return COEFFS::name(); }

  /// @brief Constructor
  ThreeSstarT(const std::string& name) : ThreeSstarImplementation(name)
  {
    m_coeffs = create_component<ThreeSstarCoeffs>("coeffs");
    m_coeffs->set( COEFFS() );
    cf3_always_assert(m_time_step_computer);
    if (m_coeffs->cfl() > 0 && m_time_step_computer->options().check("cfl") )
      m_time_step_computer->options().set("cfl",m_coeffs->cfl());
    std::cout << "constructed" << std::endl;
  }
};

/////////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // sdm
} // cf3

#endif // cf3_sdm_solver_erkls_ThreeSstar_hpp
