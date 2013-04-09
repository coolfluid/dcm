/// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test module for cf3::sdm::equations::advectiondiffusion"


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
#include "cf3/mesh/Domain.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/SimpleMeshGenerator.hpp"
#include "cf3/mesh/MeshTransformer.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/actions/LoadBalance.hpp"
#include "cf3/sdm/Model.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/PDESolver.hpp"
#include "cf3/solver/ComputeRHS.hpp"
#include "cf3/solver/BC.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/solver/TimeStepping.hpp"
#include "cf3/solver/TimeStepComputer.hpp"
#include "cf3/sdm/equations/advectiondiffusion/AdvectionDiffusion1D.hpp"

using namespace boost;
using namespace boost::assign;
using namespace cf3;
using namespace cf3::math;
using namespace cf3::common;
using namespace cf3::common::PE;
using namespace cf3::mesh;
using namespace cf3::mesh::actions;
using namespace cf3::solver;
using namespace cf3::sdm;

struct sdm_AdvectionDiffusionTests_Fixture
{
  /// common setup for each test case
  sdm_AdvectionDiffusionTests_Fixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~sdm_AdvectionDiffusionTests_Fixture()
  {
  }
  /// possibly common functions used on the tests below


  /// common values accessed by all tests goes here
  int    m_argc;
  char** m_argv;

};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( sdm_AdvectionDiffusionTests_TestSuite, sdm_AdvectionDiffusionTests_Fixture )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( init_mpi )
{
  PE::Comm::instance().init(m_argc,m_argv);
  Core::instance().environment().options().set("log_level", (Uint)INFO);
  Core::instance().environment().options().set("exception_backtrace", true);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_advection_1d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("adv1d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 1;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,100));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);
  
  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> pde = model->add_pde( /*name*/  "advectiondiffusion",
                                            /*type*/  "cf3.sdm.equations.advectiondiffusion.AdvectionDiffusion1D",
                                            /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP1" );
  pde->options().set("a",1.);
  pde->options().set("mu",0.);

  std::vector< Handle<Component> > bc_regions;
  bc_regions.push_back( mesh->access_component("topology/xneg") );
  bc_regions.push_back( mesh->access_component("topology/xpos") );
  // Handle<solver::BC> bc_mirror = pde->add_bc("mirror","cf3.sdm.equations.advectiondiffusion.BCMirror1D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  for (Uint n=0; n<pde->fields()->size(); ++n)
  {
    Real& x = pde->solution()->coordinates()[n][0];
    Real& q = pde->solution()->array()[n][0];
    if (x < 2.5)
      q = 5.;
    else
      q = -5.;
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDESolver> solver = model->add_solver( pde,
                                                        "cf3.sdm.solver.erkls.TwoSstar",
                                                        "cf3.solver.ImposeCFL" );
  solver->options().set("order",1);
  solver->time_step_computer()->options().set("cfl",1.);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",5.);
  model->time_stepping()->options().set("time_step",1.);

  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:advection1d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,pde->solution()->uri()));
  }

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  std::vector<URI> fields;
  fields.push_back(pde->solution()->uri());
  mesh->write_mesh(URI("file:advection1d.plt"), fields);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_diffusion_1d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("diff1d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 1;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,33));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);
  
  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> pde = model->add_pde( /*name*/  "advectiondiffusion",
                                            /*type*/  "cf3.sdm.equations.advectiondiffusion.AdvectionDiffusion1D",
                                            /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP2" );
  pde->options().set("a",0.);
  pde->options().set("mu",1.);

  std::vector< Handle<Component> > bc_regions;
  bc_regions.push_back( mesh->access_component("topology/xneg") );
  bc_regions.push_back( mesh->access_component("topology/xpos") );
  // Handle<solver::BC> bc_mirror = pde->add_bc("mirror","cf3.sdm.equations.advectiondiffusion.BCMirror1D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  for (Uint n=0; n<pde->fields()->size(); ++n)
  {
    Real& x = pde->solution()->coordinates()[n][0];
    Real& q = pde->solution()->array()[n][0];
    if (x < 5.)
      q = 5.;
    else
      q = -5.;
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDESolver> solver = model->add_solver( pde,
                                                        "cf3.sdm.solver.erkls.TwoSstar",
                                                        "cf3.solver.ImposeCFL" );
  solver->options().set("order",3);
  solver->time_step_computer()->options().set("cfl",0.25);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------
  
  model->time_stepping()->options().set("end_time",5.);
  model->time_stepping()->options().set("time_step",1.);
  
  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:diffusion1d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,pde->solution()->uri()));
  }

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  std::vector<URI> fields;
  fields.push_back(pde->solution()->uri());
  mesh->write_mesh(URI("file:diffusion1d.plt"), fields);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_advection_2d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("adv2d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 2;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,20));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> pde = model->add_pde( /*name*/  "advectiondiffusion",
                                            /*type*/  "cf3.sdm.equations.advectiondiffusion.AdvectionDiffusion2D",
                                            /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP0" );
  std::vector<Real> a(dim);
  a[XX] = 0.;
  a[YY] = 1.;
  pde->options().set("a",a);
  pde->options().set("mu",0.);

  std::vector< Handle<Component> > bc_regions;
//  bc_regions.push_back( mesh->access_component("topology/xneg") );
//  bc_regions.push_back( mesh->access_component("topology/xpos") );
  // Handle<solver::BC> bc_mirror = pde->add_bc("mirror","cf3.sdm.equations.advectiondiffusion.BCMirror1D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  for (Uint n=0; n<pde->fields()->size(); ++n)
  {
    Real& x = pde->solution()->coordinates()[n][XX];
    Real& y = pde->solution()->coordinates()[n][YY];
    Real& q = pde->solution()->array()[n][0];
    if (y < 2.5)
      q = 5.;
    else
      q = -5.;
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDESolver> solver = model->add_solver( pde,
                                                        "cf3.sdm.solver.erkls.TwoSstar",
                                                        "cf3.solver.ImposeCFL" );
  solver->options().set("order",1);
  solver->time_step_computer()->options().set("cfl",1.);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",5.);
  model->time_stepping()->options().set("time_step",1.);

  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:advection2d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,pde->solution()->uri()));
  }

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  std::vector<URI> fields;
  fields.push_back(pde->solution()->uri());
  mesh->write_mesh(URI("file:advection2d.plt"), fields);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_diffusion_2d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("diff2d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 2;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,20));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> pde = model->add_pde( /*name*/  "advectiondiffusion",
                                            /*type*/  "cf3.sdm.equations.advectiondiffusion.AdvectionDiffusion2D",
                                            /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP1" );
  std::vector<Real> a(dim);
  a[XX] = 0.;
  a[YY] = 0.;
  pde->options().set("a",a);
  pde->options().set("mu",1.);

  std::vector< Handle<Component> > bc_regions;
//  bc_regions.push_back( mesh->access_component("topology/xneg") );
//  bc_regions.push_back( mesh->access_component("topology/xpos") );
  // Handle<solver::BC> bc_mirror = pde->add_bc("mirror","cf3.sdm.equations.advectiondiffusion.BCMirror1D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  for (Uint n=0; n<pde->fields()->size(); ++n)
  {
    Real& x = pde->solution()->coordinates()[n][XX];
    Real& y = pde->solution()->coordinates()[n][YY];
    Real& q = pde->solution()->array()[n][0];
    if (x+y < 10.0)
      q = 5.;
    else
      q = -5.;
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDESolver> solver = model->add_solver( pde,
                                                        "cf3.sdm.solver.erkls.TwoSstar",
                                                        "cf3.solver.ImposeCFL" );
  solver->options().set("order",2);
  solver->time_step_computer()->options().set("cfl",0.5);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",5.);
  model->time_stepping()->options().set("time_step",1.);

  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:diffusion2d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,pde->solution()->uri()));
  }

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  std::vector<URI> fields;
  fields.push_back(pde->solution()->uri());
  mesh->write_mesh(URI("file:diffusion2d.plt"), fields);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
