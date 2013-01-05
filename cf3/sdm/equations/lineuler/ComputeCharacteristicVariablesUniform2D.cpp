// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/common/Builder.hpp"
#include "cf3/common/OptionComponent.hpp"
#include "cf3/common/OptionT.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/OptionArray.hpp"
#include "cf3/common/Log.hpp"
#include "cf3/common/PE/Comm.hpp"
#include "cf3/math/Consts.hpp"
#include "cf3/mesh/Field.hpp"
#include "cf3/mesh/Dictionary.hpp"
#include "cf3/mesh/FieldManager.hpp"
#include "cf3/mesh/Mesh.hpp"
#include "cf3/mesh/Region.hpp"
#include "cf3/mesh/Cells.hpp"
#include "cf3/mesh/Space.hpp"
#include "cf3/mesh/Connectivity.hpp"
#include "cf3/mesh/actions/ComputeFieldGradient.hpp"
#include "cf3/physics/lineuler/lineuler2d/Functions.hpp"
#include "cf3/sdm/equations/lineuler/ComputeCharacteristicVariablesUniform2D.hpp"


/////////////////////////////////////////////////////////////////////////////////////

using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::math::Consts;

namespace cf3 {
namespace sdm {
namespace equations {
namespace lineuler {

///////////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ComputeCharacteristicVariablesUniform2D, common::Action, LibLinEuler >
   ComputeCharacteristicVariablesUniform2D_builder;

///////////////////////////////////////////////////////////////////////////////////////

ComputeCharacteristicVariablesUniform2D::ComputeCharacteristicVariablesUniform2D ( const std::string& name ) :
  common::Action(name)
{
  mark_basic();

  std::vector<Real> normal(2);
  normal[XX]=1.;
  normal[YY]=0.;
  options().add("normal",normal).description("characteristic normal").mark_basic();
  options().add("field",m_cons)
    .link_to(&m_cons)
    .description("conservative field")
    .mark_basic();
  options().add("c0",m_c0).link_to(&m_c0).mark_basic();
}

////////////////////////////////////////////////////////////////////////////////

void ComputeCharacteristicVariablesUniform2D::execute()
{
  physics::lineuler::lineuler2d::ColVector_NDIM n;
  std::vector<Real> normal = options().value< std::vector<Real> >("normal");
  n[XX]=normal[XX];
  n[YY]=normal[YY];
  physics::lineuler::lineuler2d::RowVector_NEQS cons;
  physics::lineuler::lineuler2d::RowVector_NEQS W;
  
  if (is_null(m_char))
  {
    if (Handle<Component> found = m_cons->dict().get_child("char"))
      m_char = found->handle<Field>();
    else
      m_char = m_cons->dict().create_field("char",6).handle<Field>();
  }

  if (is_null(m_gradn_char))
  {
    if (Handle<Component> found = m_cons->dict().get_child("gradn_char"))
      m_gradn_char = found->handle<Field>();
    else
      m_gradn_char = m_cons->dict().create_field("gradn_char",12u).handle<Field>();
  }

  boost_foreach(const Handle<Entities>& elements, m_char->entities_range()  )
  {
    const Space& space = m_char->space(*elements);
    Uint nvars = m_char->row_size();
    for (Uint elem=0; elem<elements->size(); ++elem)
    {
      mesh::Connectivity::ConstRow nodes = space.connectivity()[elem];
      boost_foreach(Uint idx, nodes)
      {
        for (Uint v=0; v<nvars; ++v)
          cons[v] = m_cons->array()[idx][v];
        
        physics::lineuler::lineuler2d::cons_to_char(cons,n,m_c0,W);

        // S, Omega, Aplus, Amin, A, omega
        m_char->array()[idx][0] = W[0];
        m_char->array()[idx][1] = W[1];
        m_char->array()[idx][2] = W[2];
        m_char->array()[idx][3] = W[3];
        m_char->array()[idx][4] = 0.5*(W[2]+W[3]);
        m_char->array()[idx][5] = 0.5*(W[2]-W[3]);
      }
    }
  }
  boost::shared_ptr< mesh::actions::ComputeFieldGradient > compute_grad 
      = allocate_component< mesh::actions::ComputeFieldGradient >("compute_grad");
  compute_grad->options().set("field",m_char);
  compute_grad->options().set("field_gradient",m_gradn_char);
  compute_grad->options().set("normal",normal);
  compute_grad->execute();
}

////////////////////////////////////////////////////////////////////////////////////

} // lineuler
} // equations
} // sdm
} // cf3

////////////////////////////////////////////////////////////////////////////////////

