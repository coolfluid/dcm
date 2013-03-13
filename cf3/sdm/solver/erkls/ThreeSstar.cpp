// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/std/vector.hpp>

#include "cf3/common/Signal.hpp"
#include "cf3/common/Log.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionComponent.hpp"
#include "cf3/common/ActionDirector.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/common/Group.hpp"

#include "cf3/math/VariablesDescriptor.hpp"

#include "cf3/solver/Time.hpp"
#include "cf3/solver/TimeStepComputer.hpp"
#include "cf3/solver/Solver.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/ComputeRHS.hpp"

#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/FieldManager.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Connectivity.hpp"

#include "cf3/sdm/solver/erkls/ThreeSstar.hpp"

using namespace boost::assign;
using namespace cf3::common;
using namespace cf3::common::XML;
using namespace cf3::solver;
using namespace cf3::mesh;

namespace cf3 {
namespace sdm {
namespace solver {
namespace erkls {

///////////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ThreeSstar, common::Action, LibERKLS > ThreeSstar_Builder;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ThreeSstarImplementation::ThreeSstarImplementation ( const std::string& name ) :
  cf3::solver::PDESolver(name)
{
  options()["pde"].attach_trigger( boost::bind( &ThreeSstarImplementation::create_fields , this) );
}

///////////////////////////////////////////////////////////////////////////////////////

void ThreeSstarImplementation::create_fields()
{ 
  Handle<Component> found;

  if ( found = m_pde->fields()->get_child("backup") )
    m_backup = found->handle<Field>();
  else
    m_backup = m_pde->fields()->create_field("backup",m_pde->nb_eqs()).handle<Field>();

  if ( found = m_pde->fields()->get_child("S2") )
    m_S2 = found->handle<Field>();
  else
    m_S2 = m_pde->fields()->create_field("S2",m_pde->nb_eqs()).handle<Field>();

  if ( found = m_pde->fields()->get_child("ws") )
    m_ws = found->handle<Field>();
  else
    m_ws = m_pde->fields()->create_field("ws").handle<Field>();

  if ( found = m_pde->fields()->get_child("dt") )
    m_dt = found->handle<Field>();
  else
    m_dt = m_pde->fields()->create_field("dt").handle<Field>();
}

///////////////////////////////////////////////////////////////////////////////////////

void ThreeSstarImplementation::step()
{
  if ( is_null(m_pde) )          throw SetupError( FromHere(), "pde not configured" );
  if ( is_null(m_pde->time()) )  throw SetupError(FromHere(), "Time was not set");
  if ( is_null(m_pde->solution()) )  throw SetupError(FromHere(), "Solution was not set");
  if ( is_null(m_time_step_computer) ) throw SetupError(FromHere(), "Time step computer was not set");

  Time& time = *m_pde->time();

  m_time_step_computer->options().set("wave_speed",m_ws);
  m_time_step_computer->options().set("time_step",m_dt);
  m_pde->rhs_computer()->options().set("wave_speed",m_ws);
  m_pde->rhs_computer()->options().set("rhs",m_pde->rhs());


  Field& S1 = *m_pde->solution();;
  Field& S2 = *m_S2;
  Field& S3 = *m_backup;
  Field& R  = *m_pde->rhs();
  Field& H  = *m_dt;

  const Real T0 = time.current_time();
  const Uint nb_eqs = m_pde->nb_eqs();
  
  Real dt = 0; 
  for (Uint stage=0; stage<m_coeffs->nb_stages(); ++stage)
  {
    // Set time and iteration for this stage
    time.dt() = dt;
    time.current_time() = T0 + m_coeffs->c(stage) * dt;

    // Set boundary condition
    m_pde->bc()->execute();
    // Compute right-hand-side
    m_pde->rhs_computer()->execute();

    // Compute time step and backup solution
    if (stage == 0)
    {
      m_time_step_computer->execute();
      S2 = 0.;
      S3 = S1;
      dt = time.dt();
    }

    // Do post-processing to the rhs
    if ( is_not_null( m_pre_update ) ) m_pre_update->execute();


    /// // Use convention indexes start at 1
    /// S1 := U(t=n)   S2 := 0   S3 := U(t=n)
    /// for i = 2:m+1 do
    ///     S2 := S2 + delta(i-1)*S1
    ///     S1 := gamma(i,1)*S1 + gamma(i,2)*S2 + gamma(i,3)*S3 + beta(i,i-1)*dt*F(S1)
    /// end
    /// U(t=n+1) = S1
    /// // for error_estimate, use:
    ///     S2 := 1/sum(delta) * (S2 + delta(m+1)*S1 + delta(m+2)*S3

    for (Uint pt=0; pt<S1.size(); ++pt)
    {
      for (Uint eq=0; eq<nb_eqs; ++eq)
      {
        S2[pt][eq] =   S2[pt][eq] + m_coeffs->delta(stage)*S1[pt][eq];
        S1[pt][eq] =   m_coeffs->gamma1(stage)*S1[pt][eq]
                     + m_coeffs->gamma2(stage)*S2[pt][eq]
                     + m_coeffs->gamma3(stage)*S3[pt][eq]
                     + m_coeffs->beta(stage)*H[pt][0]*R[pt][eq];
      }
    }

    // U (=S1) has now been updated

    S1.synchronize();
    
    // Do post-processing to the stage solution
    if ( is_not_null( m_post_update ) ) m_post_update->execute();
  }
  time.current_time() = T0;
}

////////////////////////////////////////////////////////////////////////////////

ThreeSstar::ThreeSstar ( const std::string& name ) :
  ThreeSstarImplementation(name)
{
  m_coeffs = create_component<ThreeSstarCoeffs>("coeffs");
}

////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // sdm
} // cf3
