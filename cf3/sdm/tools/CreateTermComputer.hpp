// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_tools_CreateTermComputer_hpp
#define cf3_sdm_tools_CreateTermComputer_hpp

#include "cf3/sdm/tools/LibTools.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace sdm {
namespace tools {

////////////////////////////////////////////////////////////////////////////////

class sdm_tools_API CreateTermComputer : public common::Component {

public: // functions
  /// Contructor
  /// @param name of the component
  CreateTermComputer ( const std::string& name );

  /// Virtual destructor
  virtual ~CreateTermComputer() {}

  /// Get the class name
  static std::string type_name () { return "CreateTermComputer"; }

private:
  
  void signal_run(common::SignalArgs& args);
  void signature_run(common::SignalArgs& args);
};

////////////////////////////////////////////////////////////////////////////////

} // tools
} // sdm
} // cf3

#endif // cf3_sdm_tools_CreateTermComputer_hpp
