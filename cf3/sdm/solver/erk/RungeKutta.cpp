// Copyright (C) 2010-2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/bind.hpp>
#include "cf3/common/Log.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/PropertyList.hpp"
#include "cf3/common/ActionDirector.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/solver/TimeStepComputer.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/solver/ComputeRHS.hpp"
#include "cf3/sdm/solver/erk/RungeKutta.hpp"
#include "cf3/sdm/solver/erk/Types.hpp"

using namespace cf3::common;
using namespace cf3::common::XML;
using namespace cf3::solver;
using namespace cf3::mesh;

namespace cf3 {
namespace sdm {
namespace solver {
namespace erk {

////////////////////////////////////////////////////////////////////////////////

RungeKuttaImplementation::RungeKuttaImplementation ( const std::string& name ) :
  cf3::solver::PDESolver(name)
{

  options()["pde"].attach_trigger( boost::bind( &RungeKuttaImplementation::create_fields , this) );

}

///////////////////////////////////////////////////////////////////////////////////////

void RungeKuttaImplementation::create_fields()
{ 
  Handle<Component> found;

  if ( found = m_pde->fields()->get_child("backup") )
    m_backup = found->handle<Field>();
  else
    m_backup = m_pde->fields()->create_field("backup",m_pde->nb_eqs()).handle<Field>();

  Uint nb_stages = m_butcher->nb_stages();
  if (m_rhs_stages.size() != nb_stages)
    m_rhs_stages.resize(nb_stages);
  
  for (Uint i=0; i<nb_stages; ++i)
  {
    if ( found = m_pde->fields()->get_child("rhs_"+to_str(i)) )
      m_rhs_stages[i] = found->handle<Field>();
    else
      m_rhs_stages[i] = m_pde->fields()->create_field("rhs_"+to_str(i),m_pde->nb_eqs()).handle<Field>();
  }

  if ( found = m_pde->fields()->get_child("dt") )
    m_dt = found->handle<Field>();
  else
    m_dt = m_pde->fields()->create_field("dt").handle<Field>();
}

///////////////////////////////////////////////////////////////////////////////////////

void RungeKuttaImplementation::step()
{  
  const ButcherTableau& butcher = *m_butcher;
  butcher.check_throw();
  const Uint nb_stages = butcher.nb_stages();
  const Uint last_stage = nb_stages-1;

  if ( is_null(m_pde) )          throw SetupError( FromHere(), "pde not configured" );
  if ( is_null(m_pde->time()) )  throw SetupError(FromHere(), "Time was not set");
  if ( is_null(m_pde->solution()) )  throw SetupError(FromHere(), "Solution was not set");
  if ( is_null(m_time_step_computer) ) throw SetupError(FromHere(), "Time step computer was not set");

  Time& time = *m_pde->time();

  m_time_step_computer->options().set("wave_speed",m_pde->wave_speed());
  m_time_step_computer->options().set("time_step",m_dt);

  Field& U  = *m_pde->solution();
  Field& R  = *m_pde->rhs();
  Field& H  = *m_dt;
  Field& U0 = *m_backup;
  const Real T0 = time.current_time();
  const Uint nb_eqs = m_pde->nb_eqs();
  
  R = 0.;
  Real dt = 0; 
  for (Uint stage=0; stage<nb_stages; ++stage)
  {
    // Set time and iteration for this stage
    time.dt() = dt;
    time.current_time() = T0 + m_butcher->c(stage) * dt;


    // Set boundary condition
    m_pde->bc()->execute();
    // Compute right-hand-side
    m_pde->rhs_computer()->compute_rhs(*m_rhs_stages[stage],*m_pde->wave_speed());

    // Compute time step and backup solution
    if (stage == 0)
    {
      m_time_step_computer->execute();
      U0 = U;
      dt = time.dt();
    }

    // Do post-processing to the rhs
    if ( is_not_null( m_pre_update ) ) m_pre_update->execute();

    if (stage != last_stage)  // update solution for next stage
    {
      Uint next_stage = stage+1;
      /// U(s+1) = U0 + h * sum( a(s+1,j) * R(j) )
      U = U0;
      for (Uint j=0; j<next_stage; ++j)
      {
        Real a_j = butcher.a(next_stage,j);
        if (a_j != 0)
        {
          const Field& R_j = (*m_rhs_stages[j]);  // R from previous stages
          for (Uint pt=0; pt<U.size(); ++pt)
          {
            for (Uint v=0; v<U.row_size(); ++v)
            {
              U[pt][v] += H[pt][0] * a_j * R_j[pt][v];
            }
          }
        }
      }
    }
    else // weighted average of all stages forms final solution
    {
      /// U(n+1) = U(n) + h * sum( bj * Rj )
      /// R = sum( bj * Rj )
      U = U0;
      for (Uint j=0; j<nb_stages; ++j)
      {
        if (butcher.b(j)!=0)
        {
          const Field& R_j = (*m_rhs_stages[j]);  // R from previous stages
          for (Uint pt=0; pt<U.size(); ++pt)
          {
            for (Uint v=0; v<U.row_size(); ++v)
            {
              R[pt][v] += butcher.b(j) * R_j[pt][v];
              U[pt][v] += H[pt][0] * R[pt][v];
            }
          }
        }
      }
    }
    U.synchronize();
    
    // Do post-processing to the stage solution
    if ( is_not_null( m_post_update ) ) m_post_update->execute();
  }
  time.current_time() = T0;
}

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < RungeKutta, RungeKuttaImplementation, LibERK > RungeKutta_Builder;

////////////////////////////////////////////////////////////////////////////////

RungeKutta::RungeKutta ( const std::string& name ) : RungeKuttaImplementation(name)
{
  options().add("order", 4)
      .description("Order of the Runge-Kutta integration (default = RK44)\n"
                   "NOTE: This overrides any existing Butcher tableau")
      .pretty_name("RK order")
      .attach_trigger( boost::bind( &RungeKutta::config_butcher_tableau, this ) )
      .mark_basic();

  m_butcher = create_component<ButcherTableau>("butcher_tableau");
  config_butcher_tableau();
}

////////////////////////////////////////////////////////////////////////////////

// Some default coefficients that are configured with option "order"
void RungeKutta::config_butcher_tableau()
{
  const Uint order = options().value<Uint>("order");
  switch (order)
  {
    case 1: // set to forward Euler
      m_butcher->set( butcher_tableau::ForwardEuler() );
      break;
    case 2: // 2-stage 2nd-order method
      m_butcher->set( butcher_tableau::Heun2() );
      break;
    case 3: // Classic RK33
      m_butcher->set( butcher_tableau::ClassicRK33() );
      break;
    case 4: // Classic RK44
      m_butcher->set( butcher_tableau::ClassicRK44() );
      break;
    case 5: // RKF65
      m_butcher->set( butcher_tableau::RKF65() );
      break;
    default:
      CFwarn << "Cannot configure order " << order << ". Using ClassicRK44 instead." << CFendl;
      m_butcher->set( butcher_tableau::ClassicRK44() );
      break;
  }
  CFinfo << "Used Butcher tableau:\n" << m_butcher->str() << CFendl;
}

////////////////////////////////////////////////////////////////////////////////

} // erk
} // solver
} // sdm
} // cf3
