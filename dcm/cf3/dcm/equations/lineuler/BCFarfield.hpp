// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_equations_lineuler_BCFarfield_hpp
#define cf3_dcm_equations_lineuler_BCFarfield_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/core/BC.hpp"
#include "cf3/dcm/equations/lineuler/LibLinEuler.hpp"

namespace cf3 {
namespace dcm {
namespace equations {
namespace lineuler {

////////////////////////////////////////////////////////////////////////////////

class dcm_equations_lineuler_API BCFarfield2D : public dcm::core::BC<2,4>
{
public:
  BCFarfield2D(const std::string& name) : dcm::core::BC<2,4>(name) {}
  virtual ~BCFarfield2D() {}
  static std::string type_name() { return "BCFarfield2D"; }
  
  virtual void compute_boundary_solution( const RowVector_NEQS& inner_solution,
                                          const ColVector_NDIM& coords,
                                          const ColVector_NDIM& face_normal,
                                          RowVector_NEQS& boundary_solution )
  {
    boundary_solution.setZero();
  }
};

////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_lineuler_BCFarfield_hpp
