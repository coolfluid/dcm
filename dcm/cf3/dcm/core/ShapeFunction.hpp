// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_core_ShapeFunction_hpp
#define cf3_dcm_core_ShapeFunction_hpp

#include "cf3/common/BoostArray.hpp"

#include "cf3/math/MatrixTypes.hpp"

#include "cf3/mesh/ShapeFunction.hpp"
#include "cf3/mesh/GeoShape.hpp"

#include "cf3/dcm/core/LibCore.hpp"

namespace cf3 {
namespace dcm {
namespace core {

/// @brief Spectral Finite Difference shape function base class
///
/// SFD shape functions are comprised of 1D shape functions, in every direction of the
/// element dimensionality. The total shape function is then the tensorial product of these
/// 1D shape functions.
/// Therefore the only possible SFD element types are Lines (1D), Quadrilaterals(2D), Hexahedrals(3D)
class dcm_core_API ShapeFunction  : public mesh::ShapeFunction {

public:

  /// Constructor
  ShapeFunction(const std::string& name = type_name());

  /// Type name
  static std::string type_name() { return "ShapeFunction"; }

  // Concrete implementation
  virtual const RealMatrix& local_coordinates() const { return sol_pts(); }

  // Concrete implementation
  virtual const RealMatrix& mononomial_coefficients() const;

  // Concrete implementation
  virtual const RealMatrix& mononomial_exponents() const;

  // Concrete implementation
  virtual RealRowVector value(const RealVector& local_coordinate) const;

  // Concrete implementation
  virtual RealMatrix gradient(const RealVector& local_coordinate) const;

  // Concrete implementation
  virtual Uint nb_nodes() const { return nb_sol_pts(); }

  /// Compute the derivative to a given orientation using the values in flx_pts
  virtual void compute_flux_value(const Uint orientation, const RealVector& local_coordinate, RealRowVector& value) const = 0;

  /// Compute the derivative to a given orientation using the values in flx_pts
  virtual void compute_flux_derivative(const Uint orientation, const RealVector& local_coordinate, RealVector& derivative) const = 0;

  /// Number of solution points
  virtual Uint nb_sol_pts() const = 0;

  /// Number of flux points
  virtual Uint nb_flx_pts() const = 0;

  /// Solution point coordinates (rows are coordinates)
  virtual const RealMatrix& sol_pts() const = 0;

  /// Flux point coordinates (rows are coordinates)
  virtual const RealMatrix& flx_pts() const = 0;

  /// Direction flux point contributes to
  virtual Uint flx_pt_dir(const Uint flx_pt) const = 0;

  /// List of interior flux points
  virtual const std::vector<Uint>& interior_flx_pts() const = 0;

  /// List of flux points in a given face
  virtual const std::vector<Uint>& face_flx_pts(const Uint face_idx, const Uint orientation, const Uint rotation) const = 0;

  /// Sign to be multiplied with the flux computed in flx_pt
  virtual const Real& flx_pt_sign(const Uint flx_pt) const = 0;

};

////////////////////////////////////////////////////////////////////////////////

inline RealRowVector ShapeFunction::value(const RealVector& local_coordinate) const
{
  RealRowVector values(nb_sol_pts());
  compute_value(local_coordinate,values);
  return values;
}

// Concrete implementation
inline RealMatrix ShapeFunction::gradient(const RealVector& local_coordinate) const
{
  RealMatrix grad(dimensionality(),nb_sol_pts());
  compute_gradient(local_coordinate,grad);
  return grad;
}

////////////////////////////////////////////////////////////////////////////////

} // core
} // dcm
} // cf3

#endif // cf3_dcm_core_ShapeFunction_hpp
