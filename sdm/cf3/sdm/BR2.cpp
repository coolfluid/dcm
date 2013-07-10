// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"

#include "cf3/sdm/BR2.hpp"
#include "cf3/dcm/equations/advectiondiffusion/Terms.hpp"
#include "cf3/dcm/equations/euler/Terms.hpp"
#include "cf3/dcm/equations/navierstokes/Terms.hpp"
#include "cf3/dcm/equations/lineuler/Terms.hpp"
#include "cf3/common/OptionList.hpp"

using namespace cf3::dcm::equations;

namespace cf3 {
namespace sdm {

template <typename TERM>
struct BR2Builder : cf3::common::ComponentBuilder<cf3::sdm::BR2< TERM >, cf3::solver::TermComputer, cf3::sdm::LibSDM >
{
  BR2Builder( const std::string& term_computer_name ) :
    cf3::common::ComponentBuilder< cf3::sdm::BR2< TERM >, cf3::solver::TermComputer, cf3::sdm::LibSDM >( term_computer_name ) {}
};

BR2Builder<advectiondiffusion::RightHandSide1D>  ad_rhs1d("cf3.sdm.br2_advectiondiffusion_RightHandSide1D");
BR2Builder<advectiondiffusion::RightHandSide2D>  ad_rhs2d("cf3.sdm.br2_advectiondiffusion_RightHandSide2D");

BR2Builder<euler::RightHandSide1D>  eu_rhs1d("cf3.sdm.br2_euler_RightHandSide1D");
BR2Builder<euler::RightHandSide2D>  eu_rhs2d("cf3.sdm.br2_euler_RightHandSide2D");

BR2Builder<navierstokes::RightHandSide2D>  ns_rhs2d("cf3.sdm.br2_navierstokes_RightHandSide2D");
BR2Builder<navierstokes::Convection2D>     ns_conv2d ("cf3.sdm.br2_navierstokes_Convection2D");
BR2Builder<navierstokes::Diffusion2D>      ns_diff2d ("cf3.sdm.br2_navierstokes_Diffusion2D");

BR2Builder<lineuler::TermsUniform2D>           lee_rhs_uniform2d("cf3.sdm.br2_lineuler_TermsUniform2D");
BR2Builder<lineuler::SourceMonopoleUniform2D>  lee_monopole_uniform2d("cf3.sdm.br2_lineuler_SourceMonopoleUniform2D");
BR2Builder<lineuler::SourceMonopoleUniform3D>  lee_monopole_uniform3d("cf3.sdm.br2_lineuler_SourceMonopoleUniform3D");

}
}

