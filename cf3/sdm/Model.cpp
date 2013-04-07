// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/function.hpp>
#include <boost/bind.hpp>
//#include <boost/algorithm/string.hpp>
#include "cf3/common/Builder.hpp"
#include "cf3/common/Log.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/PropertyList.hpp"
//#include "common/EventHandler.hpp"
#include "cf3/common/Group.hpp"
//#include "common/Core.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/common/OptionArray.hpp"
#include "cf3/common/Signal.hpp"

#include "cf3/math/VariablesDescriptor.hpp"

#include "cf3/mesh/actions/InitFieldFunction.hpp"
#include "cf3/mesh/actions/CreateField.hpp"
#include "cf3/mesh/actions/BuildFaces.hpp"

#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/DiscontinuousDictionary.hpp"
#include "cf3/mesh/Domain.hpp"
#include "cf3/mesh/Entities.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Space.hpp"


#include "cf3/solver/History.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/solver/TimeStepping.hpp"
#include "cf3/solver/actions/ProbePostProcHistory.hpp"
#include "cf3/solver/actions/Probe.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/PDESolver.hpp"
#include "cf3/solver/ComputeRHS.hpp"

#include "cf3/sdm/Model.hpp"

#include "cf3/sdm/tools/CreateTermComputer.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::mesh::actions;
using namespace cf3::solver;
using namespace cf3::solver::actions;


namespace cf3 {
namespace sdm {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < Model, common::Component, LibSDM > Model_Builder;

////////////////////////////////////////////////////////////////////////////////

Model::Model ( const std::string& name  ) :
  common::Component ( name )
{
  mark_basic();
  // properties

  properties()["brief"] = std::string("Spectral Finite Difference Simulation model");
  properties()["description"] = std::string("Long description not available");

//  m_time = create_static_component<solver::Time>("time");
  m_time_stepping = create_static_component<solver::TimeStepping>("time_stepping");
  m_time_stepping->mark_basic();
  m_time_stepping->history()->options().set("file",URI(name+".tsv"));
  m_tools = create_static_component<Group>("tools");
  m_tools->mark_basic();
  m_domain = create_static_component<Domain>("domain");
  m_domain->mark_basic();

  regist_signal ( "add_pde" )
      .description( "Add PDE" )
      .pretty_name("Add PDE" )
      .connect   ( boost::bind ( &Model::signal_add_pde,    this, _1 ) )
      .signature ( boost::bind ( &Model::signature_add_pde, this, _1 ) );


  regist_signal ( "add_solver" )
      .description( "Add solver to solve given PDEs" )
      .pretty_name("Add Solver" )
      .connect   ( boost::bind ( &Model::signal_add_solver,    this, _1 ) )
      .signature ( boost::bind ( &Model::signature_add_solver, this, _1 ) );

  m_tools->create_component<mesh::actions::CreateField>("field_creator")->mark_basic();
  m_tools->create_component<mesh::actions::InitFieldFunction>("init_field")->mark_basic();
  m_tools->create_component<tools::CreateTermComputer>("create_term_computer")->mark_basic();
//  regist_signal("add_probe")
//      .description("Add a probe in one coordinate to inspect and log variables")
//      .pretty_name("Add Probe")
//      .connect   ( boost::bind ( &Model::signal_add_probe   , this , _1 ) )
//      .signature ( boost::bind ( &Model::signature_add_probe, this , _1 ) );

  // Set some defaults

//  m_time_stepping->post_actions() << *m_boundary_conditions;
//  (*m_pre_update) << *m_boundary_conditions;

//  m_residual_norm_computer = m_tools->create_static_component<sdm::ComputeLNorm>("compute_residual_norm");
//  m_residual_norm_computer->options().set("history",m_history);
//  (*m_post_iteration) << *m_residual_norm_computer;
}

////////////////////////////////////////////////////////////////////////////////

Model::~Model()
{
}

////////////////////////////////////////////////////////////////////////////////

Handle<mesh::Dictionary> Model::create_sd_space(const std::string& name, const Uint& order, const std::vector<Handle<Component> >& regions)
{
  if (find_components<Mesh>(*m_domain).size() == 0)
    throw SetupError(FromHere(), "Could not create solution_space because no meshes were found in "+m_domain->uri().string());

  boost_foreach(Mesh& mesh, find_components<Mesh>(*m_domain) )
  {
    build_faces(mesh);
  }

  /// @todo Support dictionary for more meshes.
  ///       This involves changing the element-finder configuration for interpolation, as it is searched for a parent of the dictionary,
  ///       to be a mesh always.
  if (find_components<Mesh>(*m_domain).size() > 1)
    throw NotImplemented(FromHere(),"Multiple mesh not supported yet");
  Handle<Mesh> mesh = find_component_ptr<Mesh>(*m_domain);
  Handle<Dictionary> dict = mesh->create_component<DiscontinuousDictionary>(name);

  std::string space_lib_name = "cf3.sdm.core.P"+to_str(order-1);
  CFinfo << "Creating Disontinuous space " << dict->uri() << " ("<<space_lib_name<<") for entities" << CFendl;
  boost_foreach(const Handle<Component>& comp, regions)
  {
    boost_foreach(const Entities& entities, find_components_recursively<Entities>( *comp ) )
    {
      CFinfo << "    -  " <<  entities.uri() << CFendl;
    }
  }

  boost_foreach(const Handle<Component>& comp, regions)
  {
    boost_foreach(Entities& entities, find_components_recursively<Entities>( *comp ) )
    {
      entities.create_space(space_lib_name+"."+entities.element_type().shape_name(),*dict);
    }
  }
  dict->build();

  boost_foreach(Mesh& mesh, find_components<Mesh>(*m_domain) )
  {
    mesh.update_structures();
  }

  return dict;
}

////////////////////////////////////////////////////////////////////////////////

void Model::build_faces(Mesh& mesh)
{
  boost::shared_ptr<BuildFaces> build_faces = allocate_component<BuildFaces>("build_inner_faces");
  build_faces->options().set("store_cell2face",true);
  build_faces->transform(mesh);
}

////////////////////////////////////////////////////////////////////////////////

void Model::signal_add_pde( common::SignalArgs& args)
{
  common::XML::SignalOptions opts(args);
  Handle<solver::PDE> pde = add_pde( opts.value<std::string>("name"),
                                     opts.value<std::string>("type"),
                                     opts.value<Uint>("order"),
                                     opts.value< std::vector< Handle<Component> > >("regions") );

  common::XML::SignalFrame reply = args.create_reply(uri());
  SignalOptions reply_options(reply);
  reply_options.add("created_component",pde->uri());
}

////////////////////////////////////////////////////////////////////////////////

void Model::signature_add_pde( common::SignalArgs& args )
{
  common::XML::SignalOptions opts(args);
  opts.add("name",std::string("pde"));
  opts.add("type",std::string("cf3.solver.PDE"));
  opts.add("order",2u);
  opts.add("regions",std::vector< Handle<Component> >(1,m_domain->handle()));
}

////////////////////////////////////////////////////////////////////////////////

void Model::signal_add_solver( common::SignalArgs& args)
{
  common::XML::SignalOptions opts(args);
  Handle<solver::PDESolver> solver = add_solver( opts.value<Handle<solver::PDE> >("pde"),
                                                 opts.value<std::string>("solver"),
                                                 opts.value<std::string>("time_step_computer") );

  common::XML::SignalFrame reply = args.create_reply(uri());
  SignalOptions reply_options(reply);
  reply_options.add("created_component",solver->uri());
}

////////////////////////////////////////////////////////////////////////////////

void Model::signature_add_solver( common::SignalArgs& args )
{
  common::XML::SignalOptions opts(args);
  opts.add("pde",Handle<solver::PDE>());
  opts.add("solver",std::string("cf3.sdm.solver.erkls.TwoSstar"));
  opts.add("time_step_computer",std::string("cf3.solver.ImposeCFL"));
}

////////////////////////////////////////////////////////////////////////////////

//void Model::signal_add_probe(common::SignalArgs& args)
//{
//  XML::SignalOptions sig_opts(args);

//  Handle<Probe> probe = m_post_iteration->create_component<Probe>(sig_opts.value<std::string>("name"));
//  probe->options().set("dict",m_solution_space);
//  probe->options().set("coordinate",sig_opts.value< std::vector<Real> >("coordinate"));
//  std::vector<std::string> functions = sig_opts.value< std::vector<std::string> >("functions");
//  boost_foreach(const std::string& function, functions)
//  {
//    std::vector<std::string> func_split;
//    boost::split(func_split,function,boost::is_any_of("="));
//    Handle<ProbePostProcessor> pp = probe->create_post_processor("func_"+func_split[0],"cf3.solver.actions.ProbePostProcFunction");
//    pp->options().set("function",function);
//  }
//  std::vector<std::string> log_vars = sig_opts.value< std::vector<std::string> >("log_variables");
//  if (log_vars.size())
//  {
//    Handle<ProbePostProcessor> pp = probe->create_post_processor("log_history","cf3.solver.actions.ProbePostProcHistory");
//    pp->options().set("history",m_history);
//    pp->options().set("variables",log_vars);
//  }
//}

//////////////////////////////////////////////////////////////////////////////////

//void Model::signature_add_probe(common::SignalArgs& args)
//{
//  XML::SignalOptions sig_opts(args);
//  sig_opts.add("name",std::string("probe"));
//  sig_opts.add("coordinate",std::vector<Real>());
//  sig_opts.add("functions",std::vector<std::string>());
//  sig_opts.add("log_variables",std::vector<std::string>());
//}

////////////////////////////////////////////////////////////////////////////////

Handle<solver::PDE> Model::add_pde(const std::string& name, const std::string& type, const Uint order)
{
  return add_pde(name,type,order,std::vector< Handle<Component> >(1,m_domain->handle()));
}

Handle<solver::PDE> Model::add_pde(const std::string& name, const std::string& type, const Uint order, const std::vector<Handle<common::Component> > &regions)
{
  Handle<Dictionary> solution_space = create_sd_space(name,order,regions);

  Handle<solver::PDE> pde = create_component(name,type)->handle<solver::PDE>();
  pde->mark_basic();
  pde->options().set("fields",solution_space);
  return pde;
}

////////////////////////////////////////////////////////////////////////////////

Handle<solver::PDESolver> Model::add_solver(const Handle<solver::PDE>& pde,
                                            const std::string& solver,
                                            const std::string& time_step_computer)
{
  Handle<solver::PDESolver> pde_solver = m_time_stepping->create_component("solve_"+pde->name(),solver)->handle<solver::PDESolver>();
  pde_solver->mark_basic();
  pde_solver->options().set("pde",pde);
  pde_solver->options().set("time_step_computer",time_step_computer);
  m_time_stepping->add_time(pde->time());
  return pde_solver;
}

////////////////////////////////////////////////////////////////////////////////

} // sdm
} // cf3
