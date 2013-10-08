// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_advectiondiffusion_BCExtrapolation_hpp
#define cf3_dcm_equations_advectiondiffusion_BCExtrapolation_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/advectiondiffusion/LibAdvectionDiffusion.hpp"
#include "cf3/common/Log.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

template <Uint NB_DIM>
class dcm_equations_advectiondiffusion_API BCExtrapolation : public dcm::core::BC<NB_DIM,1>
{
private:
  typedef typename dcm::core::BC<NB_DIM,1> BASE;
public:
  BCExtrapolation(const std::string& name) : BASE(name) {}
  virtual ~BCExtrapolation() {}
  static std::string type_name() { return "BCExtrapolation"+common::to_str(NB_DIM)+"D"; }
  
  virtual void compute_boundary_solution( const typename BASE::RowVector_NEQS& inner_solution,
                                          const typename BASE::Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const typename BASE::ColVector_NDIM& coords,
                                          const typename BASE::ColVector_NDIM& face_normal,
                                          typename BASE::RowVector_NEQS& boundary_solution,
                                          typename BASE::Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution = inner_solution;
    boundary_solution_gradient = inner_solution_gradient;
  }
};

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_advectiondiffusion_BCExtrapolation_hpp
