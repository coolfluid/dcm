/// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
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
#include "cf3/physics/euler/euler1d/Data.hpp"
#include "cf3/physics/euler/euler2d/Data.hpp"
#include "cf3/sdm/equations/euler/Euler1D.hpp"
#include "cf3/sdm/equations/euler/Euler2D.hpp"

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

struct sdm_diffusionTests_Fixture
{
  /// common setup for each test case
  sdm_diffusionTests_Fixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~sdm_diffusionTests_Fixture()
  {
  }
  /// possibly common functions used on the tests below


  /// common values accessed by all tests goes here
  int    m_argc;
  char** m_argv;

};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( sdm_diffusionTests_TestSuite, sdm_diffusionTests_Fixture )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( init_mpi )
{
  PE::Comm::instance().init(m_argc,m_argv);
  Core::instance().environment().options().set("log_level", (Uint)INFO);
  Core::instance().environment().options().set("exception_backtrace", false);
}

////////////////////////////////////////////////////////////////////////////////

#if 0
BOOST_AUTO_TEST_CASE( test_euler_1d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("model1d");

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

  Handle<solver::PDE> euler = model->add_pde( /*name*/  "euler",
                                              /*type*/  "cf3.sdm.equations.euler.Euler1D",
                                              /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP4" );
  euler->options().set("gamma",1.4);
  euler->options().set("R",287.05);
  euler->options().set("riemann_solver",std::string("Roe"));

  std::vector< Handle<Component> > bc_regions;
  bc_regions.push_back( mesh->access_component("topology/xneg") );
  bc_regions.push_back( mesh->access_component("topology/xpos") );
  Handle<solver::BC> bc_mirror = euler->add_bc("mirror","cf3.sdm.equations.euler.BCMirror1D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  physics::euler::euler1d::RowVector_NEQS prim_left, prim_right;
  physics::euler::euler1d::Data pL, pR;
  pL.gamma=1.4; pL.R=287.05;
  pR.gamma=1.4; pR.R=287.05;
  prim_left  << 4.696, 0, 404400; pL.compute_from_primitive(prim_left);
  prim_right << 1.408, 0, 101100; pR.compute_from_primitive(prim_right);

  for (Uint n=0; n<euler->fields()->size(); ++n)
  {
    if (euler->solution()->coordinates()[n][XX] < 5.)
    {
      for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
        euler->solution()->array()[n][eq] = pL.cons[eq];
    }
    else
    {
      for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
        euler->solution()->array()[n][eq] = pR.cons[eq];
    }
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

//  Handle<solver::PDESolver> solver = model->add_solver( euler,
//                                                        "cf3.sdm.solver.erkls.TwoSstar",
//                                                        "cf3.solver.ImposeCFL" );
//  solver->options().set("order",4);
//  solver->time_step_computer()->options().set("cfl",0.4);

  Handle<solver::PDESolver> solver = model->add_solver( euler,
                                                        "cf3.sdm.solver.lusgs.BDF1",
                                                        "cf3.solver.ImposeCFL" );
  solver->time_step_computer()->options().set("cfl",std::string("min(4,0.05*(i+1))"));
  std::vector<Real> Qref(euler->nb_eqs());
  Qref[0] = 1.4;
  Qref[1] = 500;
  Qref[2] = 1e5;
  solver->get_child("jacobian")->options().set("reference_solution",Qref);
  solver->options().set("max_sweeps",50);
  solver->options().set("convergence_level",1e-6);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",0.008);
  model->time_stepping()->options().set("time_step",0.001);

  try {
  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:euler1d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,euler->solution()->uri()));
  }
  } catch(...) {}

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  std::vector<URI> fields;
  fields.push_back(euler->solution()->uri());
  fields.push_back(euler->fields()->get_child("rhs")->uri());
  fields.push_back(euler->fields()->get_child("residual")->uri());
  fields.push_back(euler->fields()->get_child("dQ")->uri());
  fields.push_back(euler->fields()->get_child("convergence")->uri());
  mesh->write_mesh(URI("file:euler1d.plt"), fields);
}
#endif
////////////////////////////////////////////////////////////////////////////////

#if 1
BOOST_AUTO_TEST_CASE( test_euler_2d )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("model2d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 2;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,40));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> euler = model->add_pde( /*name*/  "euler",
                                              /*type*/  "cf3.sdm.equations.euler.Euler2D",
                                              /*sf*/    "cf3.sdm.core.LegendreGaussLobattoP0" );

  euler->options().set("gamma",1.4);
  euler->options().set("R",287.05);
  euler->options().set("riemann_solver",std::string("Roe"));

  std::vector< Handle<Component> > bc_regions;
  bc_regions.push_back( mesh->access_component("topology/left") );
  bc_regions.push_back( mesh->access_component("topology/bottom") );
  bc_regions.push_back( mesh->access_component("topology/right") );
  bc_regions.push_back( mesh->access_component("topology/top") );

  Handle<solver::BC> bc_mirror = euler->add_bc("mirror","cf3.sdm.equations.euler.BCMirror2D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  physics::euler::euler2d::RowVector_NEQS prim_left, prim_right;
  physics::euler::euler2d::Data pL, pR;
  pL.gamma=1.4; pL.R=287.05;
  pR.gamma=1.4; pR.R=287.05;
  prim_left  << 4.696, 0, 0, 404400; pL.compute_from_primitive(prim_left);
  prim_right << 1.408, 0, 0, 101100; pR.compute_from_primitive(prim_right);

  boost_foreach( const Cells& cells, find_components_recursively<Cells>(*mesh))
  {
    const Space& space = euler->fields()->space(cells);

    for (Uint elem=0; elem<cells.size(); ++elem)
    {
      boost_foreach( const Uint n, space.connectivity()[elem] )
      {
        if (euler->solution()->coordinates()[n][XX] < 5. &&
            euler->solution()->coordinates()[n][YY] < 5. )
        {
          for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
            euler->solution()->array()[n][eq] = pL.cons[eq];
        }
        else
        {
          for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
            euler->solution()->array()[n][eq] = pR.cons[eq];
        }
      }
    }
  }

//  bc_extrapolation->execute();


//  for (Uint n=0; n<euler->fields()->size(); ++n)
//  {
//    if (euler->solution()->coordinates()[n][XX] < 5. &&
//        euler->solution()->coordinates()[n][YY] < 5. )
//    {
//      for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
//        euler->solution()->array()[n][eq] = pL.cons[eq];
//    }
//    else
//    {
//      for (Uint eq=0; eq<euler->nb_eqs(); ++eq)
//        euler->solution()->array()[n][eq] = pR.cons[eq];
//    }
//  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

//  Handle<solver::PDESolver> solver = model->add_solver( euler,
//                                                        "cf3.sdm.solver.optim_erkls.ERK_5_3",
//                                                        "cf3.solver.ImposeCFL" );
//  Handle<solver::PDESolver> solver = model->add_solver( euler,
//                                                        "cf3.sdm.solver.erkls.TwoSstar",
//                                                        "cf3.solver.ImposeCFL" );
 Handle<solver::PDESolver> solver = model->add_solver( euler,
                                                       "cf3.sdm.solver.erk.MidPoint",
                                                       "cf3.solver.ImposeCFL" );
 solver->time_step_computer()->options().set("cfl",0.9);
 
  // Handle<solver::PDESolver> solver = model->add_solver( euler,
  //                                                       "cf3.sdm.solver.lusgs.BDF1",
  //                                                       "cf3.solver.ImposeCFL" );
//  solver->options().set("order",4);
  // solver->time_step_computer()->options().set("cfl",std::string("min(1.,0.1*(i+1))"));
  // std::vector<Real> Qref(euler->nb_eqs());
  // Qref[0] = 1.4;
  // Qref[1] = 500;
  // Qref[2] = 500;
  // Qref[3] = 1e5;
  // solver->get_child("jacobian")->options().set("reference_solution",Qref);
  // solver->options().set("max_sweeps",50);
  // solver->options().set("convergence_level",1e-5);
  // solver->options().set("recompute_lhs_frequency",1);
  // solver->options().set("debug",true);

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",0.008);
  model->time_stepping()->options().set("time_step",0.001);

  while ( model->time_stepping()->not_finished() )
  {
    model->time_stepping()->do_step();
    mesh->write_mesh(URI("file:euler2d_"+model->time_stepping()->options()["step"].value_str()+".plt"), std::vector<URI>(1,euler->solution()->uri()));
  }

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------
  CFinfo << "Writing final solution" << CFendl;
  std::vector<URI> fields;
  fields.push_back(euler->solution()->uri());
  // fields.push_back(euler->fields()->get_child_checked("rhs")->uri());
  // fields.push_back(euler->fields()->get_child_checked("residual")->uri());
  mesh->write_mesh(URI("file:euler2d.plt"), fields);
  mesh->write_mesh(URI("file:euler2d.msh"), fields);
  CFinfo << mesh->tree() << CFendl;
}
#endif
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
