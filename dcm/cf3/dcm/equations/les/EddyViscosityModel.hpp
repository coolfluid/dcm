// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// @file EddyViscosityModel.hpp
/// @brief Eddy viscosity models for LES
///
/// This file contains the EddyViscosityModel base class and implementations:
///   - Smagorinsky Model
///   - Wall-adapting local eddy-viscosity (WALE) model
///
/// @author Willem Deconinck


#ifndef cf3_dcm_equations_les_EddyViscosityModel_hpp
#define cf3_dcm_equations_les_EddyViscosityModel_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/dcm/equations/les/LibLES.hpp"
#include "cf3/physics/MatrixTypes.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace dcm {
namespace equations {
namespace les {

////////////////////////////////////////////////////////////////////////////////

/// @brief Base class for LES subfilter-scale Eddy Viscosity Model
///
/// This component computes the modeled values for
///   - sub-filter-scale kinetic energy
///   - sub-filter-scale kinematic viscosity a.k.a. eddy viscosity
///   - sub-filter-scale heat conduction
///
/// @author Willem Deconinck

class dcm_equations_les_API EddyViscosityModel : public common::Component
{
public: 

  /// @brief Constructor
  EddyViscosityModel( const std::string& name ) : common::Component(name) {}
  
  /// @brief Destructor
  virtual ~EddyViscosityModel() {}
  
  static std::string type_name() { return "EddyViscosityModel"; }

public: // functions

  /// Compute model in 2d
  virtual void compute( const Real& rho, const RealVector2& U, const RealMatrix2& grad_U,
                        const Real& filter_width );

  /// Compute model in 3d
  virtual void compute( const Real& rho, const RealVector3& U, const RealMatrix3& grad_U,
                        const Real& filter_width );

  /// @brief sub filter scale kinematic viscosity
  const Real& sfs_viscosity() const { return nuT; }

  /// @brief sub filter scale kinetic energy
  const Real& sfs_kinetic_energy() const { return k_sfs; }

  /// @brief sub filter scale heat conduction coefficient
  const Real& sfs_heat_conduction() const { return kappaT; }
  
protected: // data    
  
  Real nuT;    //!<  Turbulent eddy viscosity
  Real k_sfs;  //!<  Sub-filter-scale kinetic energy
  Real kappaT; //!<  Turbulent heat conduction
};

////////////////////////////////////////////////////////////////////////////////

/// @brief Smagorinsky Model
/// @author Willem Deconinck
class dcm_equations_les_API Smagorinsky : public EddyViscosityModel
{
public: 
  
  /// @brief Constructor
  Smagorinsky( const std::string& name );
  
  /// @brief Destructor
  virtual ~Smagorinsky() {}
  
  static std::string type_name() { return "Smagorinsky"; }

public: // functions

  // 2D implementation
  virtual void compute( const Real& rho, const RealVector2& U, const RealMatrix2& grad_U,
                        const Real& filter_width );

private: // data    
  
  Real PrT;
  Real Cs;
  Real Cv;
};

////////////////////////////////////////////////////////////////////////////////

/// @brief Wall-adapting local eddy-viscosity (WALE) model
/// @author Willem Deconinck
class dcm_equations_les_API WALE : public EddyViscosityModel
{
public:

  /// @brief Constructor
  WALE( const std::string& name );

  /// @brief Destructor
  virtual ~WALE() {}

  static std::string type_name() { return "WALE"; }

public: // functions

  virtual void compute( const Real& rho, const RealVector2& U, const RealMatrix2& grad_U,
                        const Real& filter_width );

private: // data

  Real PrT;
  Real Cw;
  Real Cv;
};


} // les
} // equations
} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_equations_les_EddyViscosityModel_hpp
