// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/sdm/br2/BR2.hpp"
#include "cf3/dcm/equations/advectiondiffusion/Terms.hpp"
#include "cf3/dcm/equations/euler/Terms.hpp"
#include "cf3/dcm/equations/navierstokes/Terms.hpp"
#include "cf3/dcm/equations/lineuler/Terms.hpp"
#include "cf3/dcm/equations/les/Terms.hpp"
#include "cf3/common/OptionList.hpp"

using namespace cf3::dcm::equations;

namespace cf3 {
namespace sdm {
namespace br2 {

template <typename TERM>
struct BR2Builder : cf3::common::ComponentBuilder<BR2< TERM >, cf3::solver::TermComputer, LibBR2 >
{
  BR2Builder( const std::string& term_computer_name ) :
    cf3::common::ComponentBuilder< BR2< TERM >, cf3::solver::TermComputer, LibBR2 >( term_computer_name ) {}
};

BR2Builder<advectiondiffusion::RightHandSide1D>  ad_rhs1d("cf3.sdm.br2.advectiondiffusion_RightHandSide1D");
BR2Builder<advectiondiffusion::RightHandSide2D>  ad_rhs2d("cf3.sdm.br2.advectiondiffusion_RightHandSide2D");

BR2Builder<euler::RightHandSide1D>  eu_rhs1d("cf3.sdm.br2.euler_RightHandSide1D");
BR2Builder<euler::RightHandSide2D>  eu_rhs2d("cf3.sdm.br2.euler_RightHandSide2D");

BR2Builder<navierstokes::RightHandSide2D>  ns_rhs2d("cf3.sdm.br2.navierstokes_RightHandSide2D");
BR2Builder<navierstokes::Convection2D>     ns_conv2d ("cf3.sdm.br2.navierstokes_Convection2D");
BR2Builder<navierstokes::Diffusion2D>      ns_diff2d ("cf3.sdm.br2.navierstokes_Diffusion2D");

BR2Builder<lineuler::RightHandSide2D>          lee_rhs2d("cf3.sdm.br2.lineuler_RightHandSide2D");
BR2Builder<lineuler::SourceMonopoleUniform2D>  lee_monopole_uniform2d("cf3.sdm.br2.lineuler_SourceMonopoleUniform2D");
BR2Builder<lineuler::SourceMonopoleUniform3D>  lee_monopole_uniform3d("cf3.sdm.br2.lineuler_SourceMonopoleUniform3D");
BR2Builder<lineuler::SourceDipole2D>           lee_dipole2d("cf3.sdm.br2.lineuler_SourceDipole2D");
BR2Builder<lineuler::SourceQuadrupole2D>       lee_quadrupole2d("cf3.sdm.br2.lineuler_SourceQuadrupole2D");

BR2Builder<les::RightHandSide2D>  les_rhs2d("cf3.sdm.br2.les_RightHandSide2D");

} // br2
} // sdm
} // cf3

