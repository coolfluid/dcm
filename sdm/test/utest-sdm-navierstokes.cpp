/// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test module for cf3::dcm"


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
#include "cf3/physics/euler/euler1d/Data.hpp"
#include "cf3/physics/euler/euler2d/Data.hpp"
#include "cf3/physics/navierstokes/navierstokes2d/Data.hpp"
#include "cf3/dcm/equations/euler/Euler1D.hpp"
#include "cf3/dcm/equations/euler/Euler2D.hpp"

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
  Core::instance().environment().options().set("log_level", (Uint)INFO);
  Core::instance().environment().options().set("exception_backtrace", false);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_navierstokes_2d )
{
  // Create simulation model
  Handle<dcm::Model> model = Core::instance().root().create_component<dcm::Model>("model2d");

  // ---------------------------------------------------------------------------------------
  //      CREATE MESH
  // ---------------------------------------------------------------------------------------

  Uint dim = 2;
  Handle<Mesh>          mesh           = model->domain()->create_component<Mesh>("mesh");
  Handle<MeshGenerator> mesh_generator = model->tools()->create_component("mesh_generator","cf3.mesh.SimpleMeshGenerator")->handle<MeshGenerator>();
  mesh_generator->options().set("mesh",mesh->uri());
  mesh_generator->options().set("nb_cells",std::vector<Uint>(dim,1));
  mesh_generator->options().set("lengths",std::vector<Real>(dim,10));
  mesh_generator->execute();
  allocate_component<LoadBalance>("repartitioner")->transform(mesh);

  // ---------------------------------------------------------------------------------------
  //      CREATE PHYSICS
  // ---------------------------------------------------------------------------------------

  Handle<solver::PDE> pde = model->add_pde( /*name*/  "navierstokes",
                                            /*type*/  "cf3.dcm.equations.navierstokes.NavierStokes2D",
                                            /*sf*/    "cf3.dcm.core.LegendreGaussEndP5" );
  pde->options().set("gamma",2.);
  pde->options().set("R",1.);
  pde->options().set("riemann_solver",std::string("Roe"));

  pde->add_term("rhs","cf3.sdm.br2_navierstokes_RightHandSide2D");

  std::vector< Handle<Component> > bc_regions;
  bc_regions.push_back( mesh->access_component("topology/left") );
  bc_regions.push_back( mesh->access_component("topology/bottom") );
  bc_regions.push_back( mesh->access_component("topology/right") );
  bc_regions.push_back( mesh->access_component("topology/top") );

  Handle<solver::BC> bc_wall = pde->add_bc("mirror","cf3.dcm.equations.navierstokes.BCWall2D",bc_regions);

  // ---------------------------------------------------------------------------------------
  //      INITIALISE SOLUTION
  // ---------------------------------------------------------------------------------------

  const Real gamma = 2.;
  const Real R = 1.;
  const Real dTdx = 2.;
  const Real dTdy = 3.;
  const Real u = 1.;
  const Real v = 2.;
  const Real T0 = 273.;

  boost_foreach( const Cells& cells, find_components_recursively<Cells>(*mesh))
  {
    const Space& space = pde->fields()->space(cells);

    for (Uint elem=0; elem<cells.size(); ++elem)
    {
      boost_foreach( const Uint n, space.connectivity()[elem] )
      {
        const Real x = pde->solution()->coordinates()[n][XX];
        const Real y = pde->solution()->coordinates()[n][YY];
        const Real rho = 1.;

        const Real T = T0+dTdx*x+dTdy*y;
        const Real p = rho * R * T;
        const Real rhoE = p/(gamma-1) + 0.5*rho*(u*u+v*v);

        pde->solution()->array()[n][0] = rho;
        pde->solution()->array()[n][1] = rho*u;
        pde->solution()->array()[n][2] = rho*v;
        pde->solution()->array()[n][3] = rhoE;
      }
    }
  }


  bc_wall->execute();

  for (Uint i=0; i<pde->bdry_solution()->size(); ++i)
  {
    const Real u_x = 0.;
    const Real u_y = 0.;
    const Real v_x = 0.;
    const Real v_y = 0.;
    const Real rho_x = 0.;
    const Real rho_y = 0.;

    const Real rho  = pde->bdry_solution()->array()[i][0];
    const Real u    = pde->bdry_solution()->array()[i][1]/rho;
    const Real v    = pde->bdry_solution()->array()[i][2]/rho;
    const Real E    = pde->bdry_solution()->array()[i][3]/rho;
    const Real p    = (gamma-1)*(rho*E - 0.5*rho*(u*u+v*v) );
    const Real T    = p/(rho*R);
    const Real rhoE_x = 1/(gamma-1)*(rho*R*dTdx+R*T*rho_x) + rho*(u*u_x+v*v_x) + 0.5*rho_x*(u*u+v*v);
    const Real rhoE_y = 1/(gamma-1)*(rho*R*dTdy+R*T*rho_y) + rho*(u*u_y+v*v_y) + 0.5*rho_y*(u*u+v*v);

    CFinfo << pde->bdry_solution_gradient()->array()[i][3] << "     " << pde->bdry_solution_gradient()->array()[i][7] << CFendl;
  }

  pde->rhs_computer()->execute();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
  PE::Comm::instance().finalize();
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
