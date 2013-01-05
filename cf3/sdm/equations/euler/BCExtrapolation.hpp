// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_equations_euler_BCExtrapolation_hpp
#define cf3_sdm_equations_euler_BCExtrapolation_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/sdm/core/BC.hpp"
#include "cf3/sdm/equations/euler/LibEuler.hpp"

namespace cf3 {
namespace sdm {
namespace equations {
namespace euler {

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API BCExtrapolation1D : public sdm::core::BC<1,3>
{
public:
  BCExtrapolation1D(const std::string& name) : sdm::core::BC<1,3>(name) {}
  virtual ~BCExtrapolation1D() {}
  static std::string type_name() { return "BCExtrapolation1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal, 
                                          RowVector_NEQS& boundary_solution )
  {
    boundary_solution = inner_solution;
  }
};

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API BCExtrapolation2D : public sdm::core::BC<2,4>
{
public:
  BCExtrapolation2D(const std::string& name) : sdm::core::BC<2,4>(name) {}
  virtual ~BCExtrapolation2D() {}
  static std::string type_name() { return "BCExtrapolation2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
  {
    boundary_solution = inner_solution;
  }
};

////////////////////////////////////////////////////////////////////////////////

class sdm_equations_euler_API BCExtrapolation3D : public sdm::core::BC<3,5>
{
public:
  BCExtrapolation3D(const std::string& name) : sdm::core::BC<3,5>(name) {}
  virtual ~BCExtrapolation3D() {}
  static std::string type_name() { return "BCExtrapolation3D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
  {
    boundary_solution = inner_solution;
  }
};

////////////////////////////////////////////////////////////////////////////////

} // euler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_sdm_equations_euler_BCExtrapolation_hpp
