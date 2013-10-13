// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file dcm/explicit/RungeKutta.hpp
/// @author Willem Deconinck
/// @author Matteo Parsani
///
/// This file includes the RungeKuttaImplementation component class,
/// defining the Explicit Runge Kutta scheme,
/// as well as the RungeKutta component class, and a templated
/// RungeKuttaT<BUTCHERTABLEAU> component class , deriving from
/// RungeKuttaImplementation.
///
/// The RungeKutta component is configurable
/// with a "order" option which will automatically configure some
/// default Butcher tableau's.
///
/// The RungeKuttaT<BUTCHERTABLEAU> component is a compile-time
/// configured RungeKuttaImplementation component class, by a templated
/// BUTCHERTABLEAU struct. These structs are defined in
/// the file dcm/explicit/Types.hpp

#ifndef cf3_dcm_solver_erk_RungeKutta_hpp
#define cf3_dcm_solver_erk_RungeKutta_hpp

#include "cf3/solver/PDESolver.hpp"
#include "cf3/dcm/solver/erk/LibERK.hpp"
#include "cf3/dcm/solver/erk/ButcherTableau.hpp"

namespace cf3 {
  namespace solver {
    class TimeStepComputer;
  }
}

namespace cf3 {
namespace dcm {
namespace solver {
namespace erk {

////////////////////////////////////////////////////////////////////////////////

/// @brief Runge-Kutta integration method
///
/// Standard Explicit Runge-Kutta integration using Butcher tableau coefficients
/// @verbatim
/// 0  |
/// c2 | a21
/// c3 | a31  a32
///  : |  :       `-.
/// cs | as1  as2  ..  as,s-1
/// ------------------------------
///    | b1   b2   ..   bs-1    bs
/// @endverbatim
/// @code
/// t0 = t
/// U0 = U
/// for i = 1 : s do
///     t := t0 + c(i)*dt
///     U := U0;
///     for j = 1:i-1 do
///         U := U + dt * a(i,j) * R(j);
///     end
///     R(i) := F( U );
/// end
/// U = U0
/// for i = 1 : s do
///     U := U + dt * b(i) * R(i)
/// end
/// t = t0
/// @endcode
///
/// Configuration of the coefficients a is a complete matrix, even though the matrix is lower-triangular
/// a:array[real]=a11,a12,...,a1s, a21,a22,...,a2s, ... , as1,as2,...,ass
/// b:array[real]=b1,b2,...,bs
/// c:array[real]=c1,c2,...,cs
///
/// @author Willem Deconinck
class dcm_solver_erk_API RungeKuttaImplementation : public cf3::solver::PDESolver {

public: // functions

  /// Contructor
  /// @param name of the component
  RungeKuttaImplementation ( const std::string& name );

  /// Virtual destructor
  virtual ~RungeKuttaImplementation() {}

  /// Get the class name
  static std::string type_name () { return "RungeKuttaImplementation"; }

  /// execute the action
  virtual void step ();

private: // functions

  virtual void create_fields();

protected:
  Handle<ButcherTableau> m_butcher;

private: // data

  Handle<mesh::Field> m_dt;
  
  // Registers necessary for general runge kutta algorithm
  Handle<mesh::Field> m_backup;              ///< U0
  std::vector< Handle<mesh::Field> > m_rhs_stages;  ///< R(i)
};

////////////////////////////////////////////////////////////////////////////////

/// @brief Explicit Runge Kutta integration method
///
/// Configuring the option "order" automatically preconfigures some default
/// Butcher tableaux:
/// - order = 1 --> forward Euler (ForwardEuler)
/// - order = 2 --> 2-stage 2nd-order Heun method (Heun2)
/// - order = 3 --> classic 3-stage 3rd-order Runge-Kutta method (ClassicRK33)
/// - order = 4 --> classic 4-stage 4th-order Runge-Kutta method (ClassicRK44)
/// - order = 5 --> 6-stage 5th-order Runge-Kutta-Fehlberg method (from the Fehlberg pair) (RKF65)
/// If nothing is configured, ClassicRK44 is assumed.
///
/// @author Willem Deconinck
class dcm_solver_erk_API RungeKutta : public RungeKuttaImplementation {

public: // functions

  /// @brief Contructor
  /// @param name of the component
  RungeKutta ( const std::string& name );

  /// @brief Get the class name
  static std::string type_name () { return "RungeKutta"; }

private: // functions

  /// @brief Configure the butcher tableau when the option "order" is configured
  virtual void config_butcher_tableau();

};

////////////////////////////////////////////////////////////////////////////////

/// @brief Runge Kutta time integration method, templated by a butcher tableau
///
/// The typename of the component will be:
/// "cf3.dcm.explicit.<name>"  with <name> to be replaced by
/// the name of the Butcher tableau. For a list of Butcher tableaux, check
/// dcm/explicit/Types.hpp
///
/// @author Willem Deconinck
template <typename BUTCHER_TABLEAU>
class RungeKuttaT : public RungeKuttaImplementation {
public:

  /// @brief Type name
  static std::string type_name () { return BUTCHER_TABLEAU::name(); }

  /// @brief Constructor
  RungeKuttaT(const std::string& name) : RungeKuttaImplementation(name)
  {
    m_butcher = create_component<ButcherTableau>("butcher_tableau");
    m_butcher->set( BUTCHER_TABLEAU() );
  }
};

////////////////////////////////////////////////////////////////////////////////

} // erk
} // solver
} // dcm
} // cf3

#endif // cf3_dcm_solver_erk_RungeKuttaImplementation_hpp
