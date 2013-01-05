// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_tools_ComputeLNorm_hpp
#define cf3_sdm_tools_ComputeLNorm_hpp

#include "cf3/common/Action.hpp"
#include "cf3/sdm/tools/LibTools.hpp"

/////////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace mesh   { class Field; }
namespace solver { class History; }
namespace sdm {
namespace tools {

class sdm_tools_API ComputeLNorm : public common::Action {

public: // functions
  /// Contructor
  /// @param name of the component
  ComputeLNorm ( const std::string& name );

  /// Virtual destructor
  virtual ~ComputeLNorm() {}

  /// Get the class name
  static std::string type_name () { return "ComputeLNorm"; }

  /// execute the action
  virtual void execute ();

  std::vector<Real> compute_norm( mesh::Field& field) const;

private:

  Uint compute_nb_rows(const mesh::Field& field) const;

  Handle<mesh::Field> m_field;

  Handle<solver::History> m_history;
};

////////////////////////////////////////////////////////////////////////////////

} // tools
} // sdm
} // cf3

#endif // cf3_sdm_tools_ComputeLNorm_hpp
