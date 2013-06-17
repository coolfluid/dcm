// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Reconstructions.hpp"
#include "cf3/solver/TermComputer.hpp"
#include "cf3/sdm/core/CombinedTermComputer.hpp"
#include "cf3/sdm/equations/lineuler/TermsUniform2D.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

common::ComponentBuilder<TermsUniform2D,solver::Term,LibLinEuler> TermsUniform2D_builder;
common::ComponentBuilder<core::CombinedTermComputer<TermsUniform2D>,solver::TermComputer,LibLinEuler> TermsUniform2DComputer_builder;

////////////////////////////////////////////////////////////////////////////////

TermsUniform2D::TermsUniform2D( const std::string& name ) :
  solver::TermBase<2,4,4,0>(name)
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .mark_basic()
      .description("Heat capacity ratio (Cp/Cv)");
  options().add("rho0",m_rho0).link_to(&m_rho0)
      .mark_basic()
      .description("Constant mean density");
  options().add("U0",m_U0).link_to(&m_U0)
      .mark_basic()
      .description("Constant mean velocity");
  options().add("p0",m_p0).link_to(&m_p0)
      .mark_basic()
      .description("Constant mean pressure");
}

void TermsUniform2D::get_variables( const mesh::Space& space,
                                    const Uint elem_idx,
                                    const ColVector_NDIM& coords,
                                    const mesh::ReconstructPoint& interpolation,
                                    const std::vector<mesh::ReconstructPoint>& gradient,
                                    const Matrix_NDIMxNDIM& jacobian,
                                    const Matrix_NDIMxNDIM& jacobian_inverse,
                                    const Real& jacobian_determinant,
                                    RowVector_NVAR& vars,
                                    RowVector_NGRAD& gradvars,
                                    Matrix_NDIMxNGRAD& gradvars_grad )
{
  mesh::Connectivity::ConstRow nodes = space.connectivity()[elem_idx];
  vars.setZero();
  boost_foreach( Uint n, interpolation.used_points() )
  {
    const Real L = interpolation.coeff(n);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += L * solution()->array()[nodes[n]][eq];
    }
  }
}

void TermsUniform2D::get_bdry_variables( const mesh::Space& space,
                                         const Uint elem_idx,
                                         const ColVector_NDIM& coords,
                                         const mesh::ReconstructPoint& interpolation,
                                         const std::vector<mesh::ReconstructPoint>& gradient,
                                         const Matrix_NDIMxNDIM& jacobian,
                                         const Matrix_NDIMxNDIM& jacobian_inverse,
                                         const Real& jacobian_determinant,
                                         RowVector_NVAR& vars,
                                         RowVector_NGRAD& gradvars,
                                         Matrix_NDIMxNGRAD& gradvars_grad )
{
  mesh::Connectivity::ConstRow nodes = space.connectivity()[elem_idx];
  vars.setZero();
  boost_foreach( Uint n, interpolation.used_points() )
  {
    const Real L = interpolation.coeff(n);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += L * bdry_solution()->array()[nodes[n]][eq];
    }
  }
}


void TermsUniform2D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.rho0=m_rho0;
  phys_data.U0[XX]=m_U0[XX];
  phys_data.U0[YY]=m_U0[YY];
  phys_data.p0=m_p0;
  phys_data.c0=std::sqrt(phys_data.gamma*phys_data.p0/phys_data.rho0);
}

void TermsUniform2D::compute_phys_data( const ColVector_NDIM& coords,
                                const RowVector_NVAR& vars,
                                const RowVector_NGRAD& gradvars,
                                const Matrix_NDIMxNGRAD& gradvars_grad,
                                DATA& phys_data )
{
  phys_data.compute_from_conservative(vars);
}

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3
