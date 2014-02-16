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
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"
#include "cf3/dcm/equations/lineuler/SourceDipole.hpp"

//////////////////////////////////////////////////////////////////////////////

using namespace cf3::math::Consts;
  
namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

//////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder<SourceDipole2D,solver::Term,LibLinEuler> SourceDipole2D_builder;

/////////////////////////////////////////////////////////////////////////////

SourceDipole2D::SourceDipole2D( const std::string& name )
  : solver::TermBase<2,4,0,0> (name),
    m_theta(0.),
    m_freq(1./30.),
    m_source_loc(2,0.),
    m_alpha(2.0),
    m_eps(0.5)
{
  options().add("freq",m_freq)
      .link_to(&m_freq)
      .mark_basic()
      .description("Frequency");

  options().add("location",m_source_loc)
      .link_to(&m_source_loc)
      .mark_basic()
      .description("Source location");
  
  options().add("width",m_alpha)
      .link_to(&m_alpha)
      .mark_basic()
      .description("Source width");
    
  options().add("amplitude",m_eps)
      .link_to(&m_eps)
      .mark_basic()
      .description("Source amplitude");

  options().add("angle",m_theta)
      .link_to(&m_theta)
      .mark_basic()
      .description("Angle of directivity in degrees");

}

/////////////////////////////////////////////////////////////////////////////

void SourceDipole2D::compute_phys_data( const ColVector_NDIM& coords,
                                               const RowVector_NVAR& vars,
                                               const RowVector_NGRAD& gradvars,
                                               const Matrix_NDIMxNGRAD& gradvars_grad,
                                               DATA& phys_data )
{
  phys_data.coords=coords;
}

/////////////////////////////////////////////////////////////////////////////

void SourceDipole2D::compute_source( const DATA& p, RowVector_NEQS& source )
{
  // angular frequency
  const Real omega = 2.*pi()*m_freq;

  // time
  const Real& t = m_time->current_time();

  // dipole velocity
  const Real u = fx(p.coords) * sin(omega*t);

  source[0] = 0.;
  source[1] = cos(m_theta*pi()/180.)*u;
  source[2] = sin(m_theta*pi()/180.)*u;
  source[3] = 0.;
}

/////////////////////////////////////////////////////////////////////////////

Real SourceDipole2D::fx(const ColVector_NDIM& coords)
{
  using namespace std;
  using namespace math::Consts;

  const Real dx = coords[XX]-m_source_loc[XX];
  const Real dy = coords[YY]-m_source_loc[YY];

  if (abs(dx) > m_alpha)
    return 0.;
  else
    return m_eps*cos( pi()/2. * dx/m_alpha ) * exp( -log(2.)/m_alpha * pow(dy,2) );
}

/////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3
