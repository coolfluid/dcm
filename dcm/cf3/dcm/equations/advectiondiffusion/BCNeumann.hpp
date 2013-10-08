// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_advectiondiffusion_BCNeumann_hpp
#define cf3_dcm_equations_advectiondiffusion_BCNeumann_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/advectiondiffusion/LibAdvectionDiffusion.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace advectiondiffusion {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCNeumann1D : public dcm::core::BC<1,1>
{
public:
  BCNeumann1D(const std::string& name) : dcm::core::BC<1,1>(name), m_dQdn(0.)
  {
    options().add("dQdn",m_dQdn).link_to(&m_dQdn).mark_basic();
  }
  virtual ~BCNeumann1D() {}
  static std::string type_name() { return "BCNeumann1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution[0] = inner_solution[0];
    boundary_solution_gradient = m_dQdn*face_normal;
  }
private:
  Real m_dQdn;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCNeumann2D : public dcm::core::BC<2,1>
{
public:
  BCNeumann2D(const std::string& name) : dcm::core::BC<2,1>(name), m_dQdn(0.)
  {
    options().add("dQdn",m_dQdn).link_to(&m_dQdn).mark_basic();
  }
  virtual ~BCNeumann2D() {}
  static std::string type_name() { return "BCNeumann2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution = inner_solution;
    const Real inner_gradient_normal = inner_solution_gradient.transpose()*face_normal;
    boundary_solution_gradient = inner_solution_gradient - (inner_gradient_normal+m_dQdn)*face_normal;
  }
private:
  Real m_dQdn;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_advectiondiffusion_API BCNeumann3D : public dcm::core::BC<3,1>
{
public:
  BCNeumann3D(const std::string& name) : dcm::core::BC<3,1>(name), m_dQdn(0.)
  {
    options().add("dQdn",m_dQdn).link_to(&m_dQdn).mark_basic();
  }
  virtual ~BCNeumann3D() {}
  static std::string type_name() { return "BCNeumann3D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const Matrix_NDIMxNEQS& inner_solution_gradient,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution,
                                          Matrix_NDIMxNEQS& boundary_solution_gradient )
  {
    boundary_solution = inner_solution;
    const Real inner_gradient_normal = inner_solution_gradient.transpose()*face_normal;
    boundary_solution_gradient = inner_solution_gradient - (inner_gradient_normal+m_dQdn)*face_normal;
  }
private:
  Real m_dQdn;
};

////////////////////////////////////////////////////////////////////////////////

} // advectiondiffusion
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_advectiondiffusion_BCNeumann_hpp
