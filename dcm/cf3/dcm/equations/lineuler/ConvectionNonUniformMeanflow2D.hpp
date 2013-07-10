// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_ConvectionNonUniformMeanflow2D_hpp
#define cf3_dcm_equations_lineuler_ConvectionNonUniformMeanflow2D_hpp

////////////////////////////////////////////////////////////////////////////////

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "cf3/dcm/ConvectiveTerm.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"
#include "Physics/LinEuler/Cons2D.hpp"



////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

struct NonUniformPhysData2D : PhysDataBase<4,2u>
{
  RealVectorNDIM U0;
  Real rho0;
  Real p0;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API ConvectionNonUniformMeanflow2D : public ConvectiveTerm< NonUniformPhysData2D >
{
private:
  typedef physics::LinEuler::Cons2D PHYS;

public:
  static std::string type_name() { return "ConvectionNonUniformMeanflow2D"; }
  ConvectionNonUniformMeanflow2D(const std::string& name) : ConvectiveTerm< PhysData >(name)
  {
    p.gamma = 1.4;
    options().add("gamma",p.gamma).mark_basic()
        .description("Specific heat reatio")
        .attach_trigger( boost::bind( &ConvectionNonUniformMeanflow2D::config_constants, this) );

    options().add("rho0",m_rho0).mark_basic()
        .description("Field containing mean density")
        .link_to(&m_rho0);

    options().add("U0",m_U0).mark_basic()
        .description("Field containing mean velocity")
        .link_to(&m_U0);


    options().add("p0",m_p0).mark_basic()
        .description("Field containing mean pressure")
        .link_to(&m_p0);


    p.rho0 = 1.;
    p.u0.setZero();
    p.P0 = 1.;
    p.c = 1.4;
    p.inv_c = 1./p.c;

    flx_pt_rho0.resize(1);
    flx_pt_U0.resize(NDIM);
    flx_pt_p0.resize(1);
  }

  virtual void compute_flx_pt_phys_data(const SFDElement& elem, const Uint flx_pt, PhysData& phys_data )
  {
    ConvectiveTerm<PhysData>::compute_flx_pt_phys_data(elem,flx_pt,phys_data);
    mesh::Field::View sol_pt_rho0 = m_rho0->view(elem.space->connectivity()[elem.idx]);
    mesh::Field::View sol_pt_U0   = m_U0->view(elem.space->connectivity()[elem.idx]);
    mesh::Field::View sol_pt_p0   = m_p0->view(elem.space->connectivity()[elem.idx]);
    //DEBUG  elem.reconstruct_from_solution_space_to_flux_points[flx_pt](sol_pt_rho0,phys_data.rho0) ;
    //DEBUG  elem.reconstruct_from_solution_space_to_flux_points[flx_pt](sol_pt_U0,phys_data.U0);
    //DEBUG  elem.reconstruct_from_solution_space_to_flux_points[flx_pt](sol_pt_p0,phys_data.p0);
  }


  virtual void compute_sol_pt_phys_data(const SFDElement& elem, const Uint sol_pt, PhysData& phys_data )
  {
    ConvectiveTerm<PhysData>::compute_sol_pt_phys_data(elem,sol_pt,phys_data);
    mesh::Field::View sol_pt_rho0 = m_rho0->view(elem.space->connectivity()[elem.idx]);
    mesh::Field::View sol_pt_U0   = m_U0->view(elem.space->connectivity()[elem.idx]);
    mesh::Field::View sol_pt_p0   = m_p0->view(elem.space->connectivity()[elem.idx]);

    phys_data.rho0   = sol_pt_rho0[sol_pt][0];
    phys_data.U0[XX] = sol_pt_U0[sol_pt][XX];
    phys_data.U0[YY] = sol_pt_U0[sol_pt][YY];
    phys_data.p0     = sol_pt_p0[sol_pt][0];
  }

  void config_constants()
  {
    p.gamma = options().option("gamma").value<Real>();
//    p.rho0  = options().option("rho0").value<Real>();
//    p.P0  = options().option("p0").value<Real>();

//    p.inv_rho0 = 1./p.rho0;

//    p.c=sqrt(p.gamma*p.P0*p.inv_rho0);
//    p.inv_c = 1./p.c;

//    std::vector<Real> U0 = options().option("U0").value<std::vector<Real> >();
//    for (Uint d=0; d<U0.size(); ++d)
//      p.u0[d] = U0[d];
  }

  virtual ~ConvectionNonUniformMeanflow2D() {}


  void set_meanflow_properties(const PhysData& phys_data)
  {
    p.rho0 = phys_data.rho0;
    p.u0 = phys_data.U0;
    p.P0 = phys_data.p0;
    p.inv_rho0 = 1./p.rho0;
    p.c=sqrt(p.gamma*p.P0*p.inv_rho0);
    p.inv_c = 1./p.c;
  }

  virtual void compute_analytical_flux(PhysData& data, const RealVectorNDIM& unit_normal,
                                       RealVectorNEQS& flux, Real& wave_speed)
  {
    set_meanflow_properties(data);
    PHYS::compute_properties(data.coord, data.solution , dummy_grads, p);
    PHYS::flux(p, unit_normal, flux);
    PHYS::flux_jacobian_eigen_values(p, unit_normal, eigenvalues);
    wave_speed = eigenvalues.cwiseAbs().maxCoeff();
  }

  virtual void compute_numerical_flux(PhysData& left, PhysData& right, const RealVectorNDIM& unit_normal,
                                      RealVectorNEQS& flux, Real& wave_speed)
  {
    set_meanflow_properties(left);

    // Compute left and right fluxes
    PHYS::compute_properties(left.coord, left.solution, dummy_grads, p);
    PHYS::flux(p , unit_normal, flux_left);

    PHYS::compute_properties(left.coord, right.solution, dummy_grads, p);
    PHYS::flux(p , unit_normal, flux_right);

    // Compute the averaged properties
    sol_avg.noalias() = 0.5*(left.solution+right.solution);
    PHYS::compute_properties(left.coord, sol_avg, dummy_grads, p);

    // Compute absolute jacobian using averaged properties
    PHYS::flux_jacobian_eigen_structure(p,unit_normal,right_eigenvectors,left_eigenvectors,eigenvalues);
    abs_jacobian.noalias() = right_eigenvectors * eigenvalues.cwiseAbs().asDiagonal() * left_eigenvectors;

    // flux = central flux - upwind flux
    flux.noalias() = 0.5*(flux_left + flux_right);
    flux.noalias() -= 0.5*abs_jacobian*(right.solution-left.solution);
    wave_speed = eigenvalues.cwiseAbs().maxCoeff();
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:

  Handle<mesh::Field> m_rho0;
  Handle<mesh::Field> m_U0;
  Handle<mesh::Field> m_p0;

  PHYS::MODEL::Properties p;
  PHYS::MODEL::Properties p_left;
  PHYS::MODEL::Properties p_right;

  PHYS::MODEL::SolM dummy_grads;
  PHYS::MODEL::GeoV dummy_coords;

  PHYS::MODEL::SolV sol_avg;

  PHYS::MODEL::SolV flux_left;
  PHYS::MODEL::SolV flux_right;

  PHYS::MODEL::SolV eigenvalues;
  PHYS::MODEL::JacM right_eigenvectors;
  PHYS::MODEL::JacM left_eigenvectors;
  PHYS::MODEL::JacM  abs_jacobian;

  std::vector<Real> flx_pt_rho0;
  std::vector<Real> flx_pt_U0;
  std::vector<Real> flx_pt_p0;
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_ConvectionNonUniformMeanflow2D_hpp
