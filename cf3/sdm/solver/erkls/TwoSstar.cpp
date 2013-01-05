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

#include "cf3/sdm/solver/erkls/TwoSstar.hpp"

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

common::ComponentBuilder < TwoSstar, common::Action, LibERKLS > TwoSstar_Builder;

///////////////////////////////////////////////////////////////////////////////////////

TwoSstarImplementation::TwoSstarImplementation ( const std::string& name ) :
  cf3::solver::PDESolver(name)
{
  options()["pde"].attach_trigger( boost::bind( &TwoSstarImplementation::create_fields , this) );
}

///////////////////////////////////////////////////////////////////////////////////////

void TwoSstarImplementation::create_fields()
{ 
  Handle<Component> found;

  if ( found = m_pde->fields()->get_child("backup") )
    m_backup = found->handle<Field>();
  else
    m_backup = m_pde->fields()->create_field("backup",m_pde->nb_eqs()).handle<Field>();

  if ( found = m_pde->fields()->get_child("rhs") )
    m_rhs = found->handle<Field>();
  else
    m_rhs = m_pde->fields()->create_field("rhs",m_pde->nb_eqs()).handle<Field>();

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

void TwoSstarImplementation::step()
{
  if ( is_null(m_pde) )          throw SetupError( FromHere(), "pde not configured" );
  if ( is_null(m_pde->time()) )  throw SetupError(FromHere(), "Time was not set");
  if ( is_null(m_pde->solution()) )  throw SetupError(FromHere(), "Solution was not set");
  if ( is_null(m_time_step_computer) ) throw SetupError(FromHere(), "Time step computer was not set");

  Time& time = *m_pde->time();

  m_time_step_computer->options().set("wave_speed",m_ws);
  m_time_step_computer->options().set("time_step",m_dt);
  m_pde->rhs()->options().set("wave_speed",m_ws);
  m_pde->rhs()->options().set("rhs",m_rhs);

  Field& U  = *m_pde->solution();
  Field& R  = *m_rhs;
  Field& H  = *m_dt;
  Field& U0 = *m_backup;
  const Real T0 = time.current_time();
  const Uint nb_eqs = m_pde->nb_eqs();
  
  Real dt = 0; 
  for (Uint stage=0; stage<m_coeffs->nb_stages(); ++stage)
  {
    // Set time and iteration for this stage
    time.dt() = dt;
    time.current_time() = T0 + m_coeffs->gamma(stage) * dt;


    // Set boundary condition
    m_pde->bc()->execute();
    // Compute right-hand-side
    m_pde->rhs()->execute();

    // Compute time step and backup solution
    if (stage == 0)
    {
      m_time_step_computer->execute();
      U0 = U;
      dt = time.dt();
    }

    // Do post-processing to the rhs
    if ( is_not_null( m_pre_update ) ) m_pre_update->execute();

    // Update solution
    const Real one_minus_alpha = 1. - m_coeffs->alpha(stage);    
    for (Uint pt=0; pt<U.size(); ++pt)
    {
      for (Uint eq=0; eq<nb_eqs; ++eq)
      {
        U[pt][eq] = one_minus_alpha*U0[pt][eq] + m_coeffs->alpha(stage)*U[pt][eq] + m_coeffs->beta(stage)*H[pt][0]*R[pt][eq];
      }
    }
    U.synchronize();
    
    // Do post-processing to the stage solution
    if ( is_not_null( m_post_update ) ) m_post_update->execute();
  }
  time.current_time() = T0;
}

////////////////////////////////////////////////////////////////////////////////


TwoSstar::TwoSstar ( const std::string& name ) :
  TwoSstarImplementation(name)
{
  options().add("order", 1u)
      .description("Order of the Runge-Kutta integration")
      .pretty_name("RK order")
      .attach_trigger( boost::bind( &TwoSstar::config_order , this ) );  
  
  m_coeffs = create_component<TwoSstarCoeffs>("coeffs");
  config_order();
}

////////////////////////////////////////////////////////////////////////////////

void TwoSstar::config_order()
{
  m_coeffs->options().set( "order", options().value<Uint>("order") );

  std::vector<Real> alpha;
  std::vector<Real> beta;
  std::vector<Real> gamma;

  // Set defaults Values
  switch ( m_coeffs->order() )
  {
     case 1: // Simple Forward Euler
       alpha += 0.;
       beta  += 1.0;
       gamma += 0.0;
       break;
       
     case 2: // R-K 2
       alpha += 0.0, 0.0;
       beta  += 0.5, 1.0;
       gamma += 0.0, 0.5;
       break;
       
     case 3:  // 3rd order TVD R-K scheme
       alpha += 0.0, 1.0/4.0, 2.0/3.0;
       beta  += 1.0, 1.0/4.0, 2.0/3.0;
       gamma += 0.0, 0.5,     1.0;
       break;
       
     case 4:    // R-K 4
       alpha += 0.0,     0.0,     0.0,     0.0;
       beta  += 1.0/4.0, 1.0/3.0, 1.0/2.0, 1.0;
       gamma += 0.0,     0.5,     0.5,     1.0;
       break;
  }
  
  if (gamma[0] != 0) throw BadValue(FromHere(),"gamma[0] must be zero for consistent time marching");
  
  m_coeffs->options().set("alpha",alpha);
  m_coeffs->options().set("beta",beta);
  m_coeffs->options().set("gamma",gamma);
}

////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // sdm
} // cf3
