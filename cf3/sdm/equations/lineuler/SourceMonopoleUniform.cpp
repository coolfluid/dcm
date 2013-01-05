// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/OptionList.hpp"
#include "cf3/common/Builder.hpp"
#include "cf3/math/Consts.hpp"
#include "cf3/math/Defs.hpp"
#include "cf3/solver/Time.hpp"
#include "cf3/sdm/core/CombinedTermComputer.hpp"
#include "cf3/sdm/equations/lineuler/LibLinEuler.hpp"
#include "cf3/sdm/equations/lineuler/SourceMonopoleUniform.hpp"

//////////////////////////////////////////////////////////////////////////////

using namespace cf3::math::Consts;
  
namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<SourceMonopoleUniform2D,solver::Term,LibLinEuler> SourceMonopoleUniform2D_builder;
common::ComponentBuilder<core::CombinedTermComputer<SourceMonopoleUniform2D>,solver::TermComputer,LibLinEuler> SourceMonopoleUniform2DComputer_builder;

common::ComponentBuilder<SourceMonopoleUniform3D,solver::Term,LibLinEuler> SourceMonopoleUniform3D_builder;
common::ComponentBuilder<core::CombinedTermComputer<SourceMonopoleUniform3D>,solver::TermComputer,LibLinEuler> SourceMonopoleUniform3DComputer_builder;

/////////////////////////////////////////////////////////////////////////////

SourceMonopoleUniform2D::SourceMonopoleUniform2D( const std::string& name )
  : solver::TermBase<2,4,0,0> (name),
    m_freq(1./30.),
    m_source_loc(2,0.),
    m_alpha(2.0),
    m_eps(0.5),
    m_gamma(1.4),
    m_rho0(1.),
    m_p0(1.)
{

  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("rho0",m_rho0).link_to(&m_rho0)
      .mark_basic()
      .description("Constant mean density");
  options().add("p0",m_p0).link_to(&m_p0)
      .mark_basic()
      .description("Constant mean pressure");

  options().add("freq",m_freq)
      .link_to(&m_freq)
      .mark_basic()
      .description("Frequency");

  options().add("source_location",m_source_loc)
      .link_to(&m_source_loc)
      .mark_basic()
      .description("Source location");
  
  options().add("alpha",m_alpha)
      .link_to(&m_alpha)
      .mark_basic()
      .description("Source width");
  
  
  options().add("epsilon",m_eps)
      .link_to(&m_eps)
      .mark_basic()
      .description("Source amplitude");
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.c02 = m_gamma*m_p0/m_rho0;
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform2D::compute_phys_data( const ColVector_NDIM& coords,
                                                 const RowVector_NVAR& vars,
                                                 const RowVector_NGRAD& gradvars,
                                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                                 DATA& phys_data )
{
  phys_data.coords=coords;
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform2D::compute_source( const DATA& p, RowVector_NEQS& source )
{
  // angular frequency
  const Real omega = 2.*pi()*m_freq;

  // time
  const Real& t = m_time->current_time();

  m_source = f(p.coords) * sin(omega*t);

  source[0] = m_source;
  source[1] = 0.;
  source[2] = 0.;
  source[3] = p.c02 * m_source; 
}

/////////////////////////////////////////////////////////////////////////////

Real SourceMonopoleUniform2D::f(const ColVector_NDIM& coords)
{
  using namespace std;
  const Real dx = coords[XX]-m_source_loc[XX];
  const Real dy = coords[YY]-m_source_loc[YY];
  return m_eps * exp( -log(2.)/m_alpha*(pow(dx,2) + pow(dy,2)) );
}

/////////////////////////////////////////////////////////////////////////////


SourceMonopoleUniform3D::SourceMonopoleUniform3D( const std::string& name )
  : solver::TermBase<3,5,0,0> (name),
    m_freq(1./30.),
    m_source_loc(3,0.),
    m_alpha(2.0),
    m_eps(0.5),
    m_gamma(1.4),
    m_rho0(1.),
    m_p0(1.)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("rho0",m_rho0).link_to(&m_rho0)
      .mark_basic()
      .description("Constant mean density");
  options().add("p0",m_p0).link_to(&m_p0)
      .mark_basic()
      .description("Constant mean pressure");

  options().add("freq",m_freq)
      .link_to(&m_freq)
      .mark_basic()
      .description("Frequency");

  options().add("source_location",m_source_loc)
      .link_to(&m_source_loc)
      .mark_basic()
      .description("Source location");
  
  options().add("alpha",m_alpha)
      .link_to(&m_alpha)
      .mark_basic()
      .description("Source width");
  
  
  options().add("epsilon",m_eps)
      .link_to(&m_eps)
      .mark_basic()
      .description("Source amplitude");
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform3D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.c02 = m_gamma*m_p0/m_rho0;
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform3D::compute_phys_data( const ColVector_NDIM& coords,
                                                 const RowVector_NVAR& vars,
                                                 const RowVector_NGRAD& gradvars,
                                                 const Matrix_NDIMxNGRAD& gradvars_grad,
                                                 DATA& phys_data )
{
  phys_data.coords=coords;
}

/////////////////////////////////////////////////////////////////////////////

void SourceMonopoleUniform3D::compute_source( const DATA& p, RowVector_NEQS& source )
{
  // angular frequency
  const Real omega = 2.*pi()*m_freq;

  // time
  const Real& t = m_time->current_time();

  m_source = f(p.coords) * sin(omega*t);

  source[0] = m_source;
  source[1] = 0.;
  source[2] = 0.;
  source[3] = 0.;
  source[4] = p.c02 * m_source; 
}

/////////////////////////////////////////////////////////////////////////////

Real SourceMonopoleUniform3D::f(const ColVector_NDIM& coords)
{
  using namespace std;
  const Real dx = coords[XX]-m_source_loc[XX];
  const Real dy = coords[YY]-m_source_loc[YY];
  const Real dz = coords[ZZ]-m_source_loc[ZZ];
  return m_eps * exp( -log(2.)/m_alpha*(pow(dx,2) + pow(dy,2) + pow(dz,2)) );
}

/////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
