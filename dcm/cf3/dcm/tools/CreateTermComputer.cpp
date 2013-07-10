// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "cf3/common/Log.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/Signal.hpp"
#include "cf3/common/XML/SignalOptions.hpp"
#include "cf3/dcm/tools/CreateTermComputer.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/Term.hpp"
#include "cf3/solver/TermComputer.hpp"

namespace cf3 {
namespace dcm {
namespace tools {

cf3::common::ComponentBuilder < CreateTermComputer, common::Component, LibTools > CreateTerm_builder;

/////////////////////////////////////////////////////////////////////////////////////

CreateTermComputer::CreateTermComputer ( const std::string& name ) :
  common::Component(name)
{
  regist_signal ( "run" )
      .description( "Create a term computer" )
      .pretty_name("Run" )
      .connect   ( boost::bind ( &CreateTermComputer::signal_run,    this, _1 ) )
      .signature ( boost::bind ( &CreateTermComputer::signature_run, this, _1 ) );
}

////////////////////////////////////////////////////////////////////////////////

void CreateTermComputer::signal_run(common::SignalArgs& args)
{
  common::XML::SignalOptions options( args );
  
  Handle<solver::PDE> pde = options.value< Handle<solver::PDE> >("pde");
  if ( is_null(pde) ) throw common::SetupError( FromHere(), "pde invalid" );
  
  std::string name = options.value<std::string>("name");
  std::string type = options.value<std::string>("type");
  
  Handle<solver::TermComputer> term_computer =
    pde->create_component<solver::TermComputer>(name+"_computer",type+"Computer");

  Handle<solver::Term> term =
    term_computer->create_component<solver::Term>(name,type);
  pde->configure(term);
  term_computer->options().set("term",term->handle());
    
  Handle<mesh::Field> field     = pde->fields()->create_field(name, pde->nb_eqs()).handle<mesh::Field>();
  Handle<mesh::Field> wavespeed = pde->fields()->create_field(name+"_ws", 1).handle<mesh::Field>();
  
  term_computer->options().set("field",field);
  term_computer->options().set("term_wave_speed_field",wavespeed);
  
  common::XML::SignalFrame reply = args.create_reply(uri());
  SignalOptions reply_options(reply);
  reply_options.add("created_component", term_computer->uri());
}

void CreateTermComputer::signature_run(common::SignalArgs& args)
{
  common::XML::SignalOptions options( args );

  options.add( "name", std::string("") )
      .description( "Name to give to the term" );
  
  options.add( "type", std::string("") )
      .description( "Type of the term" );

  options.add( "pde", Handle<solver::PDE>() )
      .description( "PDE" );
}

////////////////////////////////////////////////////////////////////////////////

} // tools
} // dcm
} // cf3