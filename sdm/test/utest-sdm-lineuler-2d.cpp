// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the ElementCaches of the
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
#include "cf3/solver/Model.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/physics/PhysModel.hpp"
#include "cf3/physics/Variables.hpp"
#include "cf3/mesh/Domain.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/FieldManager.hpp"
#include "cf3/mesh/SimpleMeshGenerator.hpp"
#include "cf3/mesh/MeshTransformer.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/actions/LoadBalance.hpp"
#include "cf3/sdm/Model.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/ComputeRHS.hpp"
#include "cf3/solver/TimeStepping.hpp"
#include "cf3/sdm/core/CombinedTermComputer.hpp"
#include "cf3/sdm/Solver.hpp"
#include "cf3/sdm/TimeIntegrationStepComputer.hpp"
#include "cf3/sdm/lineuler/LinEulerUniform2D.hpp"

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
using namespace cf3::sdm::core;

struct sdm_MPITests_Fixture
{
  /// common setup for each test case
  sdm_MPITests_Fixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~sdm_MPITests_Fixture()
  {
  }
  /// possibly common functions used on the tests below


  /// common values accessed by all tests goes here
  int    m_argc;
  char** m_argv;

};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( sdm_MPITests_TestSuite, sdm_MPITests_Fixture )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( init_mpi )
{
  PE::Comm::instance().init(m_argc,m_argv);
  Core::instance().environment().options().set("log_level", (Uint)INFO);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_LinEulerUniform2D )
{
  // Create simulation model
  Handle<sdm::Model> model = Core::instance().root().create_component<sdm::Model>("model1d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 2;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,20));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,1.));
  mesh_generator->options().set("offsets",std::vector<Real>(dim,-0.5));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> lineuler = model->add_pde( /*name*/  "lineuler",
                                                 /*type*/  "cf3.sdm.lineuler.LinEulerUniform2D",
                                                 /*order*/ 3 );

  lineuler->options().set("gamma",1.4);
  lineuler->options().set("rho0",1.);
  lineuler->options().set("p0",1.);
  std::vector<Real> U0(dim,0.);
  lineuler->options().set("U0",U0);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  Real c0 = std::sqrt(1.4);

  for (Uint n=0; n<lineuler->fields()->size(); ++n)
  {
    const Real x = lineuler->solution()->coordinates()[n][XX];
    const Real y = lineuler->solution()->coordinates()[n][YY];

    lineuler->solution()->array()[n][0] = 0.001*std::exp( -( x*x + y*y )/(0.05*0.05) );
    lineuler->solution()->array()[n][1] = 0.;
    lineuler->solution()->array()[n][2] = 0.;
    lineuler->solution()->array()[n][3] = c0*c0 * 0.001*std::exp( -( x*x + y*y )/(0.05*0.05) );
  }

  // ---------------------------------------------------------------------------------------
  //      SOLVE WITH RUNGEKUTTA
  // ---------------------------------------------------------------------------------------

  Handle<sdm::Solver> solver = model->add_solver( lineuler,
                                                  "cf3.sdm.ExplicitSolver",
                                                  "cf3.sdm.ExplicitRungeKuttaLowStorage2",
                                                  "cf3.sdm.ImposeCFL" );
  solver->scheme()->options().set("nb_stages",4);
  solver->time_step()->options().set("cfl",std::string("0.2"));

  // ---------------------------------------------------------------------------------------
  //      TIME STEPPING
  // ---------------------------------------------------------------------------------------

  model->time_stepping()->options().set("end_time",0.3);
  model->time_stepping()->options().set("time_step",0.05);

  model->time_stepping()->execute();

  // ---------------------------------------------------------------------------------------
  //      WRITE SOLUTION
  // ---------------------------------------------------------------------------------------

   mesh->write_mesh(URI("file:lineuler.msh"), std::vector<URI>(1,lineuler->solution()->uri()));

}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
