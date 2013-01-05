// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include <sstream>
#include <iomanip>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/OptionT.hpp"

#include "cf3/sdm/solver/erkls/ThreeSstarCoeffs.hpp"

using namespace cf3::common;
using namespace cf3::common::XML;

namespace cf3 {
namespace sdm {
namespace solver {
namespace erkls {

///////////////////////////////////////////////////////////////////////////////////////

ThreeSstarCoeffs::ThreeSstarCoeffs(const std::string& name)
  : common::Component(name)
{
  options().add("delta", std::vector<Real>())
      .description("Coefficients delta")
      .link_to(&m_coeffs.delta)
      .attach_trigger( boost::bind( &ThreeSstarCoeffs::configure_nb_stages, this) );
  
  options().add("beta", std::vector<Real>())
      .description("Coefficients beta")
      .link_to(&m_coeffs.beta);

  options().add("gamma1", std::vector<Real>())
      .description("Coefficients gamma1")
      .link_to(&m_coeffs.gamma1);
  
  options().add("gamma2", std::vector<Real>())
      .description("Coefficients gamma2")
      .link_to(&m_coeffs.gamma2);
  
  options().add("gamma3", std::vector<Real>())
      .description("Coefficients gamma3")
      .link_to(&m_coeffs.gamma3);
  
  options().add("c", std::vector<Real>())
      .description("Coefficients c")
      .link_to(&m_coeffs.c);

  options().add("order", m_coeffs.order)
      .link_to(&m_coeffs.order);
  
  options().add("cfl", m_coeffs.cfl)
      .link_to(&m_coeffs.cfl);
}

////////////////////////////////////////////////////////////////////////////////

void ThreeSstarCoeffs::set(const Coefficients& coeffs)
{
  options().set("delta",coeffs.delta);
  options().set("beta",coeffs.beta);
  options().set("gamma1",coeffs.gamma1);
  options().set("gamma2",coeffs.gamma2);
  options().set("gamma3",coeffs.gamma3);
  options().set("delta",coeffs.delta);
  options().set("c",coeffs.c);
  options().set("order",coeffs.order);
  options().set("cfl",coeffs.cfl);
}

////////////////////////////////////////////////////////////////////////////////

bool ThreeSstarCoeffs::check_throw() const
{
  if (m_coeffs.order == 0u)
  {
    throw SetupError( FromHere(), "order of coefficients is not configured" );
    return false;
  }
  if (m_coeffs.nb_stages == 0u)
  {
    throw SetupError( FromHere(), "nb_stages of coefficients is zero, configure coefficients \"delta\"" );
    return false;
  }
  if (m_coeffs.beta.size() != m_coeffs.nb_stages)
  {
    throw SetupError( FromHere(), "mismatch between beta and nb_stages in coefficients" );
    return false;
  }
  if (m_coeffs.gamma1.size() != m_coeffs.nb_stages)
  {
    throw SetupError( FromHere(), "mismatch between gamma1 and nb_stages in coefficients" );
    return false;
  }
  if (m_coeffs.gamma2.size() != m_coeffs.nb_stages)
  {
    throw SetupError( FromHere(), "mismatch between gamma2 and nb_stages in coefficients" );
    return false;
  }
  if (m_coeffs.gamma3.size() != m_coeffs.nb_stages)
  {
    throw SetupError( FromHere(), "mismatch between gamma3 and nb_stages in coefficients" );
    return false;
  }
  
  if (m_coeffs.c.size() != m_coeffs.nb_stages)
  {
    throw SetupError( FromHere(), "mismatch between c and nb_stages in coefficients" );
    return false;
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void ThreeSstarCoeffs::configure_nb_stages()
{
  m_coeffs.nb_stages = m_coeffs.delta.size();
}

////////////////////////////////////////////////////////////////////////////////

std::string ThreeSstarCoeffs::str() const
{
  std::stringstream ss;
  for (int i=0; i<nb_stages(); ++i)
  {
    ss << std::setw(2) << i << ": " 
       << std::setw(9) << delta(i) 
       << std::setw(9) << beta(i)
       << std::setw(9) << gamma1(i)
       << std::setw(9) << gamma2(i)
       << std::setw(9) << gamma3(i);
    ss << std::endl;
  }
  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////

} // erkls
} // solver
} // sdm
} // cf3
