// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_dcm_Model_hpp
#define cf3_dcm_Model_hpp

////////////////////////////////////////////////////////////////////////////////

#include "cf3/solver/Model.hpp"
#include "cf3/dcm/LibDCM.hpp"

namespace cf3 {
namespace common { class ActionDirector; }
namespace mesh { class Dictionary; }
namespace solver { class Time; class TimeStepping; class History; class PDE; class PDESolver;class ComputeRHS; }
namespace dcm {

////////////////////////////////////////////////////////////////////////////////

/// @brief Spectral Finite Difference iterative solver
///
/// Spectral Finite Difference solver,
/// combining a forward euler time marching scheme with
/// a high-order spectral finite difference spatial scheme
/// @author Willem Deconinck
class dcm_API Model : public common::Component {

public: // functions

  /// Contructor
  /// @param name of the component
  Model ( const std::string& name );

  /// Virtual destructor
  virtual ~Model();

  /// Get the class name
  static std::string type_name () { return "Model"; }

  Handle<solver::PDE> add_pde( const std::string& name,
                               const std::string& type,
                               const std::string& shape_function );

  Handle<solver::PDE> add_pde( const std::string& name,
                               const std::string& type,
                               const std::string& shape_function,
                               const std::vector<Handle<common::Component> > &regions);

  Handle<solver::PDESolver> add_solver( const std::string &name,
                                        const Handle<solver::PDE>& pde,
                                        const std::string& type,
                                        const std::string& time_step );

  const Handle<mesh::Domain>& domain() { return m_domain; }

  const Handle<common::Group>& tools() { return m_tools; }

  const Handle<solver::TimeStepping>& time_stepping() { return m_time_stepping; }

public: // signals

  void signal_add_pde( common::SignalArgs& args);
  void signature_add_pde( common::SignalArgs& args);

  void signal_add_solver( common::SignalArgs& args);
  void signature_add_solver( common::SignalArgs& args);

  void signal_create_space( common::SignalArgs& args);
  void signature_create_space( common::SignalArgs& args);

//  void signal_add_probe(common::SignalArgs& args);
//  void signature_add_probe(common::SignalArgs& args);

public: // functions

  Handle<mesh::Dictionary> create_space(const std::string& name, const std::string& shape_function, const std::vector<Handle<common::Component> > &regions);

  Handle<mesh::Dictionary> create_bdry_space(const std::string& name, const std::string& shape_function, const std::vector<Handle<common::Component> > &regions);

private: // functions

  void config_solution();

  void build_faces(cf3::mesh::Mesh &mesh);

private: // data

  Handle<common::Group>                     m_tools;
  Handle<solver::TimeStepping>              m_time_stepping;
  Handle<mesh::Domain>                      m_domain;
  Handle<solver::History>                   m_history;
};

////////////////////////////////////////////////////////////////////////////////

} // dcm
} // cf3

////////////////////////////////////////////////////////////////////////////////

#endif // cf3_dcm_time_integration_solver_hpp
