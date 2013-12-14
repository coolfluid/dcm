// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test module for cf3::sdm"


#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>
#include "cf3/common/Log.hpp"
#include "cf3/common/Core.hpp"
#include "cf3/common/Environment.hpp"
#include "cf3/common/OSystem.hpp"
#include "cf3/common/OSystemLayer.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/List.hpp"
#include "cf3/common/Group.hpp"
#include "cf3/common/PE/Comm.hpp"
#include "cf3/common/PE/debug.hpp"
#include "cf3/math/Consts.hpp"
#include "cf3/math/VariablesDescriptor.hpp"
#include "cf3/mesh/Domain.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/FieldManager.hpp"
#include "cf3/mesh/SimpleMeshGenerator.hpp"
#include "cf3/mesh/MeshTransformer.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/actions/LoadBalance.hpp"
#include "cf3/dcm/Model.hpp"
#include "cf3/solver/ComputeRHS.hpp"
#include "cf3/sdm/br2/BR2.hpp"
#include "Term.hpp"

using namespace boost;
using namespace boost::assign;
using namespace cf3;
using namespace cf3::math;
using namespace cf3::common;
using namespace cf3::common::PE;
using namespace cf3::mesh;
using namespace cf3::mesh::actions;
using namespace cf3::solver;
using namespace cf3::sdm::br2;
using namespace cf3::dcm::core;

struct testfixture
{
  /// common setup for each test case
  testfixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~testfixture()
  {
  }
  /// possibly common functions used on the tests below


  /// common values accessed by all tests goes here
  int    m_argc;
  char** m_argv;

};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( testsuite, testfixture )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( init_mpi )
{
  PE::Comm::instance().init(m_argc,m_argv);
  Core::instance().environment().options().set("log_level", (Uint)DEBUG);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( compute_1d )
{
  static const Uint dim = 1;

  // Create simulation
  Handle<dcm::Model> model = Core::instance().root().create_component<dcm::Model>("model"+to_str(dim)+"d");

  // Create mesh
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,50));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,3.));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  model->build_faces();
  Handle<Dictionary> solution_space = model->create_space("solution_space","cf3.dcm.core.LegendreGaussEndP3", std::vector< Handle<Component> >(1,mesh->handle()));
  Handle<Dictionary> solution_space_bdry = model->create_bdry_space("solution_space_bdry","cf3.dcm.core.LegendreGaussEndP3", std::vector< Handle<Component> >(1,mesh->handle()));

  // Create and setup physics
  typedef sdm::test::Term<dim>  TERM;

  // Create Numerics to solve the physics
  shared_ptr<ComputeRHS>  compute_rhs   = allocate_component<ComputeRHS>("compute_rhs");

  typedef BR2< TERM >  TERM_COMPUTER;
  Handle< solver::TermComputer > compute_term  = compute_rhs->create_component< TERM_COMPUTER >("term_computer")->handle<solver::TermComputer>();
  compute_term->options().set("alpha",0.);
  compute_term->configure_option_recursively("fields",solution_space);
  compute_term->configure_option_recursively("bdry_fields",solution_space_bdry);

  // Solve
  mesh::Field& rhs         = solution_space->create_field("rhs",dim);
  mesh::Field& wave_speed  = solution_space->create_field("wave_speed");
  mesh::Field& term_field  = solution_space->create_field("term",dim);
  
  // BOOST_CHECK_NO_THROW (compute_rhs->compute_rhs(rhs,wave_speed));
  compute_term->compute_term(term_field,wave_speed);

  mesh->write_mesh(URI("file:term_field_"+to_str(dim)+"d.msh"), std::vector<URI>(1,term_field.uri()));
  // mesh->write_mesh(URI("file:rhs_"+to_str(dim)+"d.msh"), std::vector<URI>(1,rhs.uri()));

  std::cout << std::boolalpha;
  std::cout << "has_convection: " << TERM::ENABLE_CONVECTION << std::endl;
  std::cout << "has_diffusion:  " << TERM::ENABLE_DIFFUSION  << std::endl;
  std::cout << "has_source:     " << TERM::ENABLE_SOURCE     << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( compute_2d )
{
  static const Uint dim = 2;

  // Create simulation
  Handle<dcm::Model> model = Core::instance().root().create_component<dcm::Model>("model"+to_str(dim)+"d");

  // Create mesh
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,10));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,3.));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);
  model->build_faces();
  Handle<Dictionary> solution_space = model->create_space("solution_space","cf3.dcm.core.LegendreGaussEndP2", std::vector< Handle<Component> >(1,mesh->handle()));
  Handle<Dictionary> solution_space_bdry = model->create_bdry_space("solution_space_bdry","cf3.dcm.core.LegendreGaussEndP3", std::vector< Handle<Component> >(1,mesh->handle()));

  // Create and setup physics
  typedef sdm::test::Term<dim>  TERM;

  // Create Numerics to solve the physics
  shared_ptr<ComputeRHS>  compute_rhs   = allocate_component<ComputeRHS>("compute_rhs");

  typedef BR2< TERM >  TERM_COMPUTER;
  Handle< solver::TermComputer > compute_term  = compute_rhs->create_component< TERM_COMPUTER >("term_computer")->handle<solver::TermComputer>();
  compute_term->configure_option_recursively("fields",solution_space);
  compute_term->configure_option_recursively("bdry_fields",solution_space_bdry);

  // Solve
  mesh::Field& rhs         = solution_space->create_field("rhs",dim);
  mesh::Field& wave_speed  = solution_space->create_field("wave_speed");
  mesh::Field& term_field  = solution_space->create_field("term",dim);
  
  // BOOST_CHECK_NO_THROW (compute_rhs->compute_rhs(rhs,wave_speed));
  compute_term->compute_term(term_field,wave_speed);

  std::vector<URI> fields;
  fields.push_back(term_field.uri());
  fields.push_back(wave_speed.uri());
  mesh->write_mesh(URI("file:term_field_"+to_str(dim)+"d.msh"), fields);
  mesh->write_mesh(URI("file:rhs_"+to_str(dim)+"d.msh"), std::vector<URI>(1,rhs.uri()));
  
  std::cout << std::boolalpha;
  std::cout << "has_convection: " << TERM::ENABLE_CONVECTION << std::endl;
  std::cout << "has_diffusion:  " << TERM::ENABLE_DIFFUSION  << std::endl;
  std::cout << "has_source:     " << TERM::ENABLE_SOURCE     << std::endl;
}

////////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE( compute_3d )
{
  static const Uint dim = 3;

  // Create simulation
  Handle<dcm::Model> model = Core::instance().root().create_component<dcm::Model>("model"+to_str(dim)+"d");

  // Create mesh
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,10));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,3.));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  model->build_faces();
  Handle<Dictionary> solution_space = model->create_space("solution_space","cf3.dcm.core.LegendreGaussEndP3", std::vector< Handle<Component> >(1,mesh->handle()));
  Handle<Dictionary> solution_space_bdry = model->create_bdry_space("solution_space_bdry","cf3.dcm.core.LegendreGaussEndP3", std::vector< Handle<Component> >(1,mesh->handle()));

  // Create and setup physics
  typedef sdm::test::Term<dim>  TERM;

  // Create Numerics to solve the physics
  shared_ptr<ComputeRHS>  compute_rhs   = allocate_component<ComputeRHS>("compute_rhs");

  typedef BR2< TERM >  TERM_COMPUTER;
  Handle< solver::TermComputer > compute_term  = compute_rhs->create_component< TERM_COMPUTER >("term_computer")->handle<solver::TermComputer>();
  compute_term->configure_option_recursively("fields",solution_space);
  compute_term->configure_option_recursively("bdry_fields",solution_space_bdry);

  // Solve
  mesh::Field& rhs         = solution_space->create_field("rhs",dim);
  mesh::Field& wave_speed  = solution_space->create_field("wave_speed");
  mesh::Field& term_field  = solution_space->create_field("term",dim);
  
  // BOOST_CHECK_NO_THROW (compute_rhs->compute_rhs(rhs,wave_speed));
  compute_term->compute_term(term_field,wave_speed);

  mesh->write_mesh(URI("file:term_field_"+to_str(dim)+"d.msh"), std::vector<URI>(1,term_field.uri()));
  // mesh->write_mesh(URI("file:rhs_"+to_str(dim)+"d.msh"), std::vector<URI>(1,rhs.uri()));

  std::cout << std::boolalpha;
  std::cout << "has_convection: " << TERM::ENABLE_CONVECTION << std::endl;
  std::cout << "has_diffusion:  " << TERM::ENABLE_DIFFUSION  << std::endl;
  std::cout << "has_source:     " << TERM::ENABLE_SOURCE     << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
