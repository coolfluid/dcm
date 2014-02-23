// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/dcm/equations/les/EddyViscosityModel.hpp"
#include "cf3/math/Defs.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

cf3::common::ComponentBuilder<Smagorinsky, cf3::common::Component, LibLES > smagorinsky_builder;
cf3::common::ComponentBuilder<WALE, cf3::common::Component, LibLES > WALE_builder;

////////////////////////////////////////////////////////////////////////////////

void EddyViscosityModel::compute(const Real& rho, const RealVector2& U, const RealMatrix2& grad_U, const Real& filter_width)
{
  throw common::NotImplemented(FromHere(), "2D-implementation not implemented");
}

void EddyViscosityModel::compute(const Real& rho, const RealVector3& U, const RealMatrix3& grad_U, const Real& filter_width)
{
  throw common::NotImplemented(FromHere(), "3D-implementation not implemented");
}

////////////////////////////////////////////////////////////////////////////////

Smagorinsky::Smagorinsky( const std::string& name ) :
  EddyViscosityModel(name),
  PrT(0.9),
  Cs(0.18),
  Cv(0.094)
{
  // Subfilter scale parameters to tune
  options().add("PrT",PrT).link_to(&PrT)
      .description("Turbulent Prandtl number")
      .mark_basic();
  options().add("Cs",Cs).link_to(&Cs)
      .description("Smagorinsky constant")
      .mark_basic();
  options().add("Cv",Cv).link_to(&Cv)
      .description("Deardorff/Yoshizawa constant")
      .mark_basic();
}


// 2D implementation
void Smagorinsky::compute( const Real& rho, const RealVector2& U, const RealMatrix2& grad_U,
                           const Real& filter_width )
{
  // Strain rate tensor
  Real Sxx = grad_U(XX,0);
  Real Sxy = 0.5*(grad_U(YY,0)+grad_U(XX,1));
  Real Syy = grad_U(YY,1);

  // Absolute strain rate tensor
  Real SijSij = Sxx*Sxx + Syy*Syy + 2*Sxy*Sxy;
  Real absS = std::sqrt( 2.* SijSij);

  // Eddy viscosity computed by Smagorinsky model
  nuT = std::pow(filter_width*Cs, 2)*absS;

  // SFS thermal conductivity
  kappaT = rho*nuT/PrT;

  // Subfilterscale kinetic energy model by Yoshizawa, Deardorff
  k_sfs = std::pow( nuT/(Cv*filter_width) , 2);
}

////////////////////////////////////////////////////////////////////////////////


WALE::WALE( const std::string& name ) :
  EddyViscosityModel(name),
  PrT(0.9),
  Cw(0.325),
  Cv(0.094)
{
  // Subfilter scale parameters to tune
  options().add("PrT",PrT).link_to(&PrT)
      .description("Turbulent Prandtl number")
      .mark_basic();
  options().add("Cw",Cw).link_to(&Cw)
      .description("WALE constant")
      .mark_basic();
  options().add("Cv",Cv).link_to(&Cv)
      .description("Deardorff/Yoshizawa constant")
      .mark_basic();
}

// 2D implementation
void WALE::compute( const Real& rho, const RealVector2& U, const RealMatrix2& grad_U,
                    const Real& filter_width )
{
  enum { DXX = 0, DYY = 1 };
  const Real epsilon = 1e-6;

  Real const& dudx = grad_U(DXX,XX);
  Real const& dudy = grad_U(DYY,XX);
  Real const& dvdx = grad_U(DXX,YY);
  Real const& dvdy = grad_U(DYY,YY);

  // temporary parameters
  const Real gXX2 = dudx*dudx + dudy*dvdx;
  const Real gXY2 = dudx*dudy + dudy*dvdy;

  const Real gYX2 = dvdx*dudy + dvdy*dvdx;
  const Real gYY2 = dvdx*dudy + dvdy*dvdy;

  const Real div = (gXX2 + gYY2)*0.5;

  const Real sXXD = gXX2 - div;//0.5*(gXX2 + gXX2) - div;
  const Real sXYD = 0.5*(gXY2 + gYX2);

  const Real sYXD = 0.5*(gYX2 + gXY2);
  const Real sYYD = gYY2 - div;//0.5*(gYY2 + gYY2) - div;

  const Real sSD = sXXD*sXXD + sXYD*sXYD +
                     sYXD*sYXD + sYYD*sYYD;

  const Real sS = (dudx*dudx + dvdy*dvdy) +
                    0.5*(dudy+dvdx)*(dudy+dvdx);

  const Real sW = std::pow(sSD,1.5)/( std::pow(sS,2.5) + std::pow(sSD,1.25) + epsilon);

  // compute the eddy viscosity
  nuT = std::pow(filter_width*Cw, 2)*sW;

  // SFS thermal conductivity
  kappaT = rho*nuT/PrT;

  // Subfilterscale kinetic energy model by Yoshizawa, Deardorff
  k_sfs = std::pow( nuT/(Cv*filter_width) , 2);
}

////////////////////////////////////////////////////////////////////////////////

} // les
} // equations
} // dcm
} // cf3
