// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_advectiondiffusion_BCDirichlet_hpp
#define cf3_dcm_equations_advectiondiffusion_BCDirichlet_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/advectiondiffusion/LibAdvectionDiffusion.hpp"
#include "cf3/common/Log.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCDirichlet1D : public dcm::core::BC<1,1>
{
public:
  BCDirichlet1D(const std::string& name) : dcm::core::BC<1,1>(name), m_Q(0.)
  {
    options().add("Q",m_Q).link_to(&m_Q).mark_basic();
  }
  virtual ~BCDirichlet1D() {}
  static std::string type_name() { return "BCDirichlet1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution[0] = m_Q;
    boundary_solution_gradient = inner_solution_gradient;
  }
private:
  Real m_Q;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCDirichlet2D : public dcm::core::BC<2,1>
{
public:
  BCDirichlet2D(const std::string& name) : dcm::core::BC<2,1>(name), m_Q(0.)
  {
    options().add("Q",m_Q).link_to(&m_Q).mark_basic();
  }
  virtual ~BCDirichlet2D() {}
  static std::string type_name() { return "BCDirichlet2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution[0] = m_Q;
    boundary_solution_gradient = inner_solution_gradient;
  }
private:
  Real m_Q;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCDirichlet3D : public dcm::core::BC<3,1>
{
public:
  BCDirichlet3D(const std::string& name) : dcm::core::BC<3,1>(name), m_Q(0.)
  {
    options().add("Q",m_Q).link_to(&m_Q).mark_basic();
  }
  virtual ~BCDirichlet3D() {}
  static std::string type_name() { return "BCDirichlet3D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution[0] = m_Q;
    boundary_solution_gradient = inner_solution_gradient;
  }
private:
  Real m_Q;
};

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_advectiondiffusion_BCDirichlet_hpp
