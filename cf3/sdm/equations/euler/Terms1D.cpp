// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/bind.hpp>
#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Reconstructions.hpp"
#include "cf3/solver/TermComputer.hpp"
#include "cf3/sdm/core/CombinedTermComputer.hpp"
#include "cf3/sdm/equations/euler/Terms1D.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

common::ComponentBuilder<Terms1D,solver::Term,LibEuler> Terms1D_builder;
common::ComponentBuilder<core::CombinedTermComputer<Terms1D>,solver::TermComputer,LibEuler> Terms1DComputer_builder;

////////////////////////////////////////////////////////////////////////////////

Terms1D::Terms1D( const std::string& name ) :
  solver::TermBase<1,3,3,0>(name),
  m_riemann_solver_type("Roe")
{
  options().add("gamma",m_gamma).link_to(&m_gamma)
      .description("Heat capacity ratio (Cp/Cv)")
      .mark_basic();
  options().add("R",m_R).link_to(&m_R)
      .description("Gas constant")
      .mark_basic();
  options().add("riemann_solver",m_riemann_solver_type).link_to(&m_riemann_solver_type)
      .description("Riemann solver")
      .attach_trigger( boost::bind( &Terms1D::config_riemann_solver, this) );
  config_riemann_solver();
}

void Terms1D::get_variables( const mesh::Space& space,
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
    const Real C = interpolation.coeff(n);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += C * solution()->array()[nodes[n]][eq];
    }
  }
}

void Terms1D::get_bdry_variables( const mesh::Space& space,
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
    const Real C = interpolation.coeff(n);
    for (Uint eq=0; eq<NEQS; ++eq)
    {
      vars[eq] += C * bdry_solution()->array()[nodes[n]][eq];
    }
  }
}

void Terms1D::set_phys_data_constants( DATA& phys_data )
{
  phys_data.gamma=m_gamma;
  phys_data.R=m_R;
}

void Terms1D::compute_phys_data( const ColVector_NDIM& coords,
                                const RowVector_NVAR& vars,
                                const RowVector_NGRAD& gradvars,
                                const Matrix_NDIMxNGRAD& gradvars_grad,
                                DATA& phys_data )
{
  phys_data.compute_from_conservative(vars);
}


void Terms1D::config_riemann_solver()
{
  std::string builder_name = m_riemann_solver_type;
  if (! boost::algorithm::starts_with(builder_name,"cf3." ) )
    builder_name = "cf3.sdm.equations.euler."+builder_name+"1D";

  if (m_riemann_solver) remove_component(*m_riemann_solver);
  m_riemann_solver = create_component("riemann_solver",builder_name)->handle< solver::RiemannSolver<Terms1D> >();
}

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3
