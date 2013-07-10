// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_core_Reconstructions_hpp
#define cf3_dcm_core_Reconstructions_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/mesh/Reconstructions.hpp"
#include "cf3/dcm/core/ShapeFunction.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

struct InterpolateInPointFromFlxPts
{
  template <typename vector_type>
  static void build_coefficients( mesh::ReconstructPoint& reconstruct,
                                  const Uint direction,
                                  const vector_type& local_coord,
                                  const Handle<dcm::core::ShapeFunction const>& sf )
  {
    reconstruct.m_N.resize(sf->nb_flx_pts());
    sf->compute_flux_value(direction,local_coord,reconstruct.m_N);
    reconstruct.construct_used_points();
    cf3_always_assert(reconstruct.used_points().size());
  }
};

////////////////////////////////////////////////////////////////////////////////

struct DerivativeInPointFromFlxPts
{
  template <typename vector_type>
  static void build_coefficients( mesh::ReconstructPoint& reconstruct,
                                  const Uint derivative_to,
                                  const vector_type& local_coord,
                                  const Handle<dcm::core::ShapeFunction const>& sf )
  {
    RealVector tmp(sf->nb_flx_pts());
    sf->compute_flux_derivative(derivative_to,local_coord,tmp);
    reconstruct.m_N = tmp;
    reconstruct.construct_used_points();
    cf3_always_assert(reconstruct.used_points().size());
  }
};

////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_core_Reconstructions_hpp
