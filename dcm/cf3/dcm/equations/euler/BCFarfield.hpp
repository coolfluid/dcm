// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_euler_BCFarfield_hpp
#define cf3_dcm_equations_euler_BCFarfield_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/euler/LibEuler.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API BCFarfield1D : public dcm::core::BC<1,3>
{
public:
  BCFarfield1D(const std::string& name);
  virtual ~BCFarfield1D() {}
  static std::string type_name() { return "BCFarfield1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution );

private:
  void compute_state();
  
private:
  Real m_rho;
  Real m_u;
  Real m_p;
  Real m_gamma;
  RowVector_NEQS state;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API BCFarfield2D : public dcm::core::BC<2,4>
{
public:
  BCFarfield2D(const std::string& name);
  virtual ~BCFarfield2D() {}
  static std::string type_name() { return "BCFarfield2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution );
  
private:
  void compute_state();
  
private:
  Real m_rho;
  std::vector<Real> m_U;
  Real m_p;
  Real m_gamma;
  RowVector_NEQS state;
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_euler_API BCFarfield3D : public dcm::core::BC<3,5>
{
public:
  BCFarfield3D(const std::string& name);
  virtual ~BCFarfield3D() {}
  static std::string type_name() { return "BCFarfield3D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );

private:
  void compute_state();
  
private:
  Real m_rho;
  std::vector<Real> m_U;
  Real m_p;
  Real m_gamma;
  RowVector_NEQS state;
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_euler_BCFarfield_hpp
