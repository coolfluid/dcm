// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>

#include "common/Log.hpp"
#include "common/Builder.hpp"
#include "common/FindComponents.hpp"
#include "common/Foreach.hpp"
#include "common/StringConversion.hpp"
#include "common/PropertyList.hpp"
#include "common/OptionList.hpp"
#include "common/OptionT.hpp"

#include "math/VariablesDescriptor.hpp"
#include "math/VectorialFunction.hpp"

#include "mesh/AInterpolator.hpp"
#include "mesh/LoadMesh.hpp"
#include "mesh/Field.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Dictionary.hpp"
#include "mesh/Space.hpp"
#include "mesh/Entities.hpp"
#include "mesh/Connectivity.hpp"

#include "cf3/dcm/tools/ComputeFieldGradientBR2.hpp"

//////////////////////////////////////////////////////////////////////////////

using namespace cf3::common;
using namespace cf3::math;
using namespace cf3::mesh;

namespace cf3 {
namespace dcm {
namespace tools {

////////////////////////////////////////////////////////////////////////////////

cf3::common::ComponentBuilder < ComputeFieldGradientBR2, mesh::MeshTransformer, cf3::dcm::tools::LibTools> ComputeFieldGradientBR2_Builder;

//////////////////////////////////////////////////////////////////////////////

ComputeFieldGradientBR2::ComputeFieldGradientBR2( const std::string& name )
: MeshTransformer(name)
{
  properties()["brief"] = std::string("ComputeFieldGradientBR2 mesh");
  std::string desc;
  desc = "  Usage: ComputeFieldGradientBR2 \n\n";

  properties()["description"] = desc;

  options().add("field",m_field).link_to(&m_field)
      .description("Field to take gradient of")
      .mark_basic();

  options().add("field_gradient",m_field_gradient).link_to(&m_field_gradient)
      .description("Output: gradient of option \"field\"")
      .mark_basic();

//  options().add("normal",m_normal).link_to(&m_normal);


  options().add("alpha",-1.)
      .description("Damping coefficient in BR2 scheme for face-gradient computation\n"
                   "If negative, alpha = 1/(P+1) is used");

  options().add("BR2",true);
}

/////////////////////////////////////////////////////////////////////////////

void ComputeFieldGradientBR2::execute()
{
  // Check correct configuration
  if (is_null(m_field)) throw SetupError(FromHere(), "Option field not set in "+uri().string());
  if (is_null(m_field_gradient)) throw SetupError(FromHere(), "Option field_gradient not set in "+uri().string());


  const Uint ndim = m_field->coordinates().row_size();


  if (m_field_gradient->row_size() != ndim*m_field->row_size())
  {
    throw SetupError(FromHere(), "Field "+m_field_gradient->uri().string()+" must have row-size of "+to_str(ndim*m_field->row_size())+". Currently it is "+to_str(m_field_gradient->row_size()));
  }
  switch(ndim)
  {
    case 1u:
      compute_gradient<1u>();
      break;
    case 2u:
      compute_gradient<2u>();
      break;
    case 3u:
      compute_gradient<3u>();
      break;
  }
}

////////////////////////////////////////////////////////////////////////////////

} // tools
} // dcm
} // cf3
