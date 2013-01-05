// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "cf3/mesh/ShapeFunctionT.hpp"
#include "cf3/mesh/LagrangeP0/Point.hpp"
#include "cf3/sdm/core/LibCore.hpp"
#include "cf3/sdm/core/LagrangeLocally1D.hpp"
#include "cf3/common/Builder.hpp"

using namespace cf3::common;

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Point<0>,mesh::ShapeFunction,LibCore>
PointP0_builder(LibCore::library_namespace()+".P0.Point");

ComponentBuilder<Point<1>,mesh::ShapeFunction,LibCore>
PointP1_builder(LibCore::library_namespace()+".P1.Point");

ComponentBuilder<Point<2>,mesh::ShapeFunction,LibCore>
PointP2_builder(LibCore::library_namespace()+".P2.Point");

ComponentBuilder<Point<3>,mesh::ShapeFunction,LibCore>
PointP3_builder(LibCore::library_namespace()+".P3.Point");

ComponentBuilder<Point<4>,mesh::ShapeFunction,LibCore>
PointP4_builder(LibCore::library_namespace()+".P4.Point");

ComponentBuilder<Point<5>,mesh::ShapeFunction,LibCore>
PointP5_builder(LibCore::library_namespace()+".P5.Point");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<LineLagrange1D<0>,mesh::ShapeFunction,LibCore>
LineP0_builder(LibCore::library_namespace()+".P0.Line");

ComponentBuilder<LineLagrange1D<1>,mesh::ShapeFunction,LibCore>
LineP1_builder(LibCore::library_namespace()+".P1.Line");

ComponentBuilder<LineLagrange1D<2>,mesh::ShapeFunction,LibCore>
LineP2_builder(LibCore::library_namespace()+".P2.Line");

ComponentBuilder<LineLagrange1D<3>,mesh::ShapeFunction,LibCore>
LineP3_builder(LibCore::library_namespace()+".P3.Line");

ComponentBuilder<LineLagrange1D<4>,mesh::ShapeFunction,LibCore>
LineP4_builder(LibCore::library_namespace()+".P4.Line");

ComponentBuilder<LineLagrange1D<5>,mesh::ShapeFunction,LibCore>
LineP5_builder(LibCore::library_namespace()+".P5.Line");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<QuadLagrange1D<0>,mesh::ShapeFunction,LibCore>
QuadP0_builder(LibCore::library_namespace()+".P0.Quad");

ComponentBuilder<QuadLagrange1D<1>,mesh::ShapeFunction,LibCore>
QuadP1_builder(LibCore::library_namespace()+".P1.Quad");

ComponentBuilder<QuadLagrange1D<2>,mesh::ShapeFunction,LibCore>
QuadP2_builder(LibCore::library_namespace()+".P2.Quad");

ComponentBuilder<QuadLagrange1D<3>,mesh::ShapeFunction,LibCore>
QuadP3_builder(LibCore::library_namespace()+".P3.Quad");

ComponentBuilder<QuadLagrange1D<4>,mesh::ShapeFunction,LibCore>
QuadP4_builder(LibCore::library_namespace()+".P4.Quad");

ComponentBuilder<QuadLagrange1D<5>,mesh::ShapeFunction,LibCore>
QuadP5_builder(LibCore::library_namespace()+".P5.Quad");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<HexaLagrange1D<0>,mesh::ShapeFunction,LibCore>
HexaP0_builder(LibCore::library_namespace()+".P0.Hexa");

ComponentBuilder<HexaLagrange1D<1>,mesh::ShapeFunction,LibCore>
HexaP1_builder(LibCore::library_namespace()+".P1.Hexa");

ComponentBuilder<HexaLagrange1D<2>,mesh::ShapeFunction,LibCore>
HexaP2_builder(LibCore::library_namespace()+".P2.Hexa");

ComponentBuilder<HexaLagrange1D<3>,mesh::ShapeFunction,LibCore>
HexaP3_builder(LibCore::library_namespace()+".P3.Hexa");

ComponentBuilder<HexaLagrange1D<4>,mesh::ShapeFunction,LibCore>
HexaP4_builder(LibCore::library_namespace()+".P4.Hexa");

ComponentBuilder<HexaLagrange1D<5>,mesh::ShapeFunction,LibCore>
HexaP5_builder(LibCore::library_namespace()+".P5.Hexa");

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

