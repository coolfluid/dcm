// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/PropertyList.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/common/FindComponents.hpp"
#include "cf3/dcm/equations/navierstokes/ComputeBoundaryLayer.hpp"
#include "cf3/dcm/core/ShapeFunction.hpp"
#include "cf3/dcm/core/FaceConnectivity.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/Faces.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/ShapeFunction.hpp"
#include "cf3/dcm/core/Reconstructions.hpp"
#include "cf3/dcm/core/Metrics.hpp"

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::dcm::core;

namespace cf3 {
namespace dcm {
namespace equations {
namespace navierstokes {

////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ComputeBoundaryLayer, common::Action, LibNavierStokes > ComputeBoundaryLayer_Builder;

////////////////////////////////////////////////////////////////////////////////

ComputeBoundaryLayer::ComputeBoundaryLayer ( const std::string& name  ) :
  common::Action ( name )
{
  // properties
  properties()["brief"] = std::string("NavierStokes 2D Partial Differential Equations");
  properties()["description"] = std::string("Component that can solve the 2D NavierStokes physics right-hand-side");

  options().add("gamma",1.4)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("R",287.05)
      .mark_basic()
      .description("Gas constant");
  options().add("nu",0.)
      .mark_basic()
      .description("Kinematic viscosity");
  options().add("velocity",m_velocity)
      .description("Velocity field")
      .mark_basic()
      .link_to(&m_velocity);
  options().add("density",m_density)
      .description("Density field")
      .mark_basic()
      .link_to(&m_density);
  options().add("y0",0.)
      .description("Wall Distance")
      .mark_basic();
  options().add("wall_regions",m_wall_regions)
      .description("Region that describes the wall")
      .mark_basic()
      .link_to(&m_wall_regions);
  options().add("yplus",m_yplus)
      .description("[out] yplus Field")
      .mark_basic()
      .link_to(&m_yplus);
  options().add("tau",m_tau)
      .description("[out] Wall shear stress Field")
      .mark_basic()
      .link_to(&m_tau);
  options().add("ustar",m_ustar)
      .description("Friction velocity Field")
      .mark_basic()
      .link_to(&m_ustar);
  options().add("yplus_min",0.)
      .description("[out] minimum value of yplus Field")
      .mark_basic();
  options().add("yplus_max",0.)
      .description("[out] maximum value of yplus Field")
      .mark_basic();
}

////////////////////////////////////////////////////////////////////////////////

ComputeBoundaryLayer::~ComputeBoundaryLayer()
{
}

////////////////////////////////////////////////////////////////////////////////

void ComputeBoundaryLayer::execute()
{
  Dictionary& dictionary = m_velocity->dict();

  if( is_null(m_velocity) ) throw SetupError(FromHere(), "velocity field not setup");
  if( is_null(m_density)  ) throw SetupError(FromHere(), "density field not setup");

  const Uint dim = dictionary.options().value<Uint>("dimension");
  const Real g = options().value<Real>("gamma");
  const Real R = options().value<Real>("R");
  Field::ArrayT& U   = m_velocity->array();
  Field::ArrayT& rho = m_density->array();
  Field::ArrayT& walldist = m_walldistance->array();

  Real nu = options().value<Real>("nu"); // kinematic viscosity
  Real y0 = options().value<Real>("y0"); // walldistance

  Real yplus_min = 1.e10;
  Real yplus_max = 0.;

  Dictionary& bdry_dict = m_yplus->dict();

  // 1) Identify inner cells and connectivity with wall, and wall unit-normal

  boost_foreach( const Handle<Region>& region, m_wall_regions )
  {
    boost_foreach( Faces& faces, find_components_recursively<Faces>(*region) )
    {
      Handle<FaceConnectivity> connected( faces.get_child("face_connectivity") );
      if (is_null(connected))
        connected = faces.create_component< FaceConnectivity >("face_connectivity");
      
      const Space& bdry_space = bdry_dict.space(faces);

      const Uint nb_faces = faces.size();
      
      for (Uint f=0; f<nb_faces; ++f)
      {
        connected->compute_face(faces, f);
        cf3_always_assert( connected->is_bdry_face() );
        Cells& cells = *connected->cells()[0].comp->handle<Cells>();
        const Uint cell_idx = connected->cells()[0].idx;
        const Space& space = dictionary.space(cells);

        const dcm::core::ShapeFunction& sf = *space.shape_function().handle<dcm::core::ShapeFunction>();
        mesh::Connectivity::ConstRow nodes = space.connectivity()[cell_idx];
        
        Handle< dcm::core::Metrics<2u> > metrics( const_cast<mesh::Space*>(&space)->get_child("metrics") );
        if (is_null(metrics))
        {
          metrics = const_cast<mesh::Space*>(&space)->create_component< dcm::core::Metrics<2u> >("metrics");
          metrics->setup_for_space(space.handle<Space>());
        }
        ElementMetrics<2u>& element_metrics = *metrics->element(cell_idx);

        const std::vector<Uint>& cell_face_flx_pts = 
          sf.face_flx_pts( connected->cells_face_nb()[0],
                           connected->cells_orientation()[0],
                           connected->cells_rotation()[0] );

        std::vector<Real>        rho_wall(cell_face_flx_pts.size(), 0.);
        std::vector<RealVector2> U_wall(cell_face_flx_pts.size());
        std::vector<RealMatrix2> grad_U_wall(cell_face_flx_pts.size());
        for( Uint face_pt=0; face_pt<cell_face_flx_pts.size(); ++face_pt)
        {
          Uint cell_face_flx_pt = cell_face_flx_pts[face_pt];
          U_wall[face_pt].setZero();
          boost_foreach( const Uint sol_pt, metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).used_points() )
          {
            const int n = nodes[sol_pt];
            const Real C = metrics->interpolation_from_sol_pts_to_flx_pt(cell_face_flx_pt).coeff(sol_pt);
            rho_wall[face_pt]   += C * rho[n][0];
            U_wall[face_pt][XX] += C * U[n][XX];
            U_wall[face_pt][YY] += C * U[n][YY];
          }
          grad_U_wall[face_pt].setZero();
          for (Uint d=0; d<2; ++d)
          {
            boost_foreach( const Uint sol_pt, metrics->gradient_from_sol_pts_to_flx_pt(cell_face_flx_pt)[d].used_points() )
            {
              const Real n = nodes[sol_pt];
              const Real D = metrics->gradient_from_sol_pts_to_flx_pt(cell_face_flx_pt)[d].coeff(sol_pt);
              grad_U_wall[face_pt](d,XX) += D * U[n][XX];
              grad_U_wall[face_pt](d,YY) += D * U[n][YY];
            }
          }
          grad_U_wall[face_pt] = element_metrics.flx_pt_Jinv(cell_face_flx_pt) * grad_U_wall[face_pt];

          Real rho_w = rho_wall[face_pt];
          Real dUxdx_w = grad_U_wall[face_pt](XX,XX);
          Real dUydx_w = grad_U_wall[face_pt](XX,YY);
          Real dUxdy_w = grad_U_wall[face_pt](YY,XX);
          Real dUydy_w = grad_U_wall[face_pt](YY,YY);

          Real nx = element_metrics.flx_pt_unit_normal(cell_face_flx_pt)[XX]*sf.flx_pt_sign(cell_face_flx_pt);
          Real ny = element_metrics.flx_pt_unit_normal(cell_face_flx_pt)[YY]*sf.flx_pt_sign(cell_face_flx_pt);

          Real mu = nu/rho_w;
          // wall shear stress tau = mu dUxdy (perpendicular to wall)
          Real tau_w = mu * ( ny * (dUxdx_w*nx+dUxdy_w*ny) -nx * (dUydx_w*nx+dUydy_w*ny) );
          // friction velocity
          Real ustar = std::sqrt(std::abs(tau_w)/rho_w);
          Real yplus = ustar * y0 / nu;
          yplus_min = std::min(yplus_min, yplus);
          yplus_max = std::max(yplus_max, yplus);

          const int bdry_node = bdry_space.connectivity()[f][face_pt];

          m_yplus->array()[ bdry_node ][0] = yplus;
          m_tau->array()[ bdry_node ][0] = tau_w;
          m_ustar->array()[ bdry_node ][0] = ustar;
        }
      }
    }
    PE::Comm::instance().all_reduce(PE::min(), &yplus_min, 1, &yplus_min);
    PE::Comm::instance().all_reduce(PE::max(), &yplus_max, 1, &yplus_max);
    options().set("yplus_min",yplus_min);
    options().set("yplus_max",yplus_max);
  }
  // 2) Compute primitive variables in inner cells

  // 3) At wall: compute gradients of temperature and velocity
      
  // 4) Compute dynamic viscosity mu at wall using vars and grads
  
  // 4) Get density from vars
  
  // 5) Compute skin friction tau at wall
     /*
      if (dim == DIM_2D)
      {  
        _tau = _muWall*
             (m_unitNormal[YY]*dot_prod(*grad_Ux,m_unitNormal) -
              m_unitNormal[XX]*dot_prod(*grad_Uy,m_unitNormal));

        // friction coefficient 
        m_Cf = _tau / (0.5*m_rhoInf*m_uInf*m_uInf);

        m_frictionForces[XX] =  m_Cf*m_normal[YY];
        m_frictionForces[YY] = -m_Cf*m_normal[XX];
      }
      else if (dim == DIM_3D)
      {
        Real divU = (2.0/3.0)*((*_gradients[m_UID])[0] + (*_gradients[m_VID])[1] + (*_gradients[m_WID])[2]) ;
  
        _tau3D(XX,XX) = _muWall*(2.0*(*_gradients[m_UID])[XX] - divU) ;
        _tau3D(YY,YY) = _muWall*(2.0*(*_gradients[m_VID])[YY] - divU) ;
        _tau3D(ZZ,ZZ) = _muWall*(2.0*(*_gradients[m_WID])[ZZ] - divU) ;
  
        _tau3D(XX,YY) = _tau3D(YY,XX) = _muWall*((*_gradients[m_UID])[YY] + (*_gradients[m_VID])[XX]) ;
        _tau3D(XX,ZZ) = _tau3D(ZZ,XX) = _muWall*((*_gradients[m_UID])[ZZ] + (*_gradients[m_WID])[XX]) ;
        _tau3D(YY,ZZ) = _tau3D(ZZ,YY) = _muWall*((*_gradients[m_VID])[ZZ] + (*_gradients[m_WID])[YY]) ;
        CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getR() * m_TInf);
        m_Cf3D = _tau3D / (0.5*rhoInf*m_uInf*m_uInf);
        //Compute the viscous force coefficients on this wall face.
        m_frictionForces = m_Cf3D * m_normal;
      }
      // forces coefficients must be still divided by the wet surface
      m_frictionForces /= m_refArea;
      */
  
  // 6) Compute yplus of first cell
      /*
      _yPlus = 0.;

      DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

      // Compute the distance to the first cell: y0
      //First get the inner state
      State* innerState = m_currFace->getState(0);
      cf_assert(!innerState->isGhost());

      CFreal y0 = wallDistance[innerState->getLocalID()];

      // Compute the y+
      // y^{+} = \frac {\sqrt{\rho} * \sqrt{\tau} * y0} {\mu}

      const CFreal refSpeed = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::V];
      const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];

      _yPlus = sqrt(_rhoWall) * sqrt(fabs(_tau)) * y0 / _muWall;
      _yPlus *= sqrt(rhoRef) * refSpeed;
      */
}

////////////////////////////////////////////////////////////////////////////////

} // navierstokes
} // equations
} // dcm
} // cf3
