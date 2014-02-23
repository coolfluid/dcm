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
#include "cf3/dcm/Model.hpp"
#include "cf3/solver/PDE.hpp"
#include "cf3/solver/PDESolver.hpp"
#include "cf3/solver/ComputeRHS.hpp"
#include "cf3/solver/BC.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/solver/TimeStepping.hpp"
#include "cf3/solver/TimeStepComputer.hpp"
#include "cf3/dcm/equations/les/EddyViscosityModel.hpp"
#include "cf3/dcm/equations/les/RightHandSide2D.hpp"
#include "cf3/dcm/equations/navierstokes/RightHandSide2D.hpp"

using namespace boost;
using namespace boost::assign;
using namespace cf3;
using namespace cf3::math;
using namespace cf3::common;
using namespace cf3::common::PE;
using namespace cf3::mesh;
using namespace cf3::mesh::actions;
using namespace cf3::solver;
using namespace cf3::dcm;
using namespace cf3::dcm::equations::les;
using namespace cf3::dcm::equations;

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

BOOST_AUTO_TEST_CASE( test_les_2d )
{

  Handle<les::RightHandSide2D> les  = Core::instance().root().create_component<les::RightHandSide2D>("les_rhs");
  Handle<navierstokes::RightHandSide2D> navierstokes = Core::instance().root().create_component<navierstokes::RightHandSide2D>("ns_rhs");

  const Real tol_pct = 1.e-5;

  les::RightHandSide2D::DATA p_les;
  les->set_phys_data_constants(p_les);
  BOOST_CHECK_CLOSE( p_les.R,  287.05, tol_pct);
  BOOST_CHECK_CLOSE( p_les.gamma,  1.4, tol_pct);
  BOOST_CHECK_CLOSE( p_les.kappa,  2.601e-2, tol_pct);
  BOOST_CHECK_CLOSE( p_les.mu,  1.806e-5, tol_pct);

  navierstokes::RightHandSide2D::DATA p_ns;
  navierstokes->set_phys_data_constants(p_ns);
  BOOST_CHECK_CLOSE( p_ns.R,  287.05, tol_pct);
  BOOST_CHECK_CLOSE( p_ns.gamma,  1.4, tol_pct);
  BOOST_CHECK_CLOSE( p_ns.kappa,  2.601e-2, tol_pct);
  BOOST_CHECK_CLOSE( p_ns.mu,  1.806e-5, tol_pct);

  navierstokes::RightHandSide2D::RowVector_NEQS solution, flux_les, flux_ns;
  navierstokes::RightHandSide2D::ColVector_NDIM normal;
  Real ws_les, ws_ns;
  solution << 1., 20., 30., 1.0e5;
  normal << 1., 0.;

  p_ns.compute_from_conservative(solution);
  p_les.compute_from_conservative(solution);

  les->compute_convective_flux(p_les, normal, flux_les, ws_les );
  navierstokes->compute_convective_flux(p_ns, normal, flux_ns, ws_ns );

  BOOST_CHECK_EQUAL( flux_les, flux_ns);
  BOOST_CHECK_EQUAL( ws_les, ws_ns);

}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
