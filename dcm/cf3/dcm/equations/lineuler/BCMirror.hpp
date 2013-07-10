// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_BCMirror_hpp
#define cf3_dcm_equations_lineuler_BCMirror_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API BCMirror1D : public dcm::core::BC<1,3>
{
public:
  BCMirror1D(const std::string& name) : dcm::core::BC<1,3>(name) {}
  virtual ~BCMirror1D() {}
  static std::string type_name() { return "BCMirror1D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API BCMirror2D : public dcm::core::BC<2,4>
{
public:
  BCMirror2D(const std::string& name) : dcm::core::BC<2,4>(name) {}
  virtual ~BCMirror2D() {}
  static std::string type_name() { return "BCMirror2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
};

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API BCMirror3D : public dcm::core::BC<3,5>
{
public:
  BCMirror3D(const std::string& name) : dcm::core::BC<3,5>(name) {}
  virtual ~BCMirror3D() {}
  static std::string type_name() { return "BCMirror3D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution );
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_BCMirror_hpp
