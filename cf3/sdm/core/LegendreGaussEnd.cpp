// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


/// @file LegendreGaussEnd.cpp
///
/// This file creates and registers shape functions for tensorial
/// elements with Legendre-Gauss-Lobatto flux-points
///
/// Shape-functions are named e.g. "cf3.sdm.core.LegendreGaussEndP2"

#include "cf3/sdm/core/Tensorial.hpp"
#include "cf3/sdm/core/LegendreGaussEnd.hpp"
#include "cf3/common/Builder.hpp"

using namespace cf3::common;

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Point<LegendreGaussEnd,0>,mesh::ShapeFunction,LibCore>
PointP0_builder(LibCore::library_namespace()+".LegendreGaussEndP0.Point");

ComponentBuilder<Point<LegendreGaussEnd,1>,mesh::ShapeFunction,LibCore>
PointP1_builder(LibCore::library_namespace()+".LegendreGaussEndP1.Point");

ComponentBuilder<Point<LegendreGaussEnd,2>,mesh::ShapeFunction,LibCore>
PointP2_builder(LibCore::library_namespace()+".LegendreGaussEndP2.Point");

ComponentBuilder<Point<LegendreGaussEnd,3>,mesh::ShapeFunction,LibCore>
PointP3_builder(LibCore::library_namespace()+".LegendreGaussEndP3.Point");

ComponentBuilder<Point<LegendreGaussEnd,4>,mesh::ShapeFunction,LibCore>
PointP4_builder(LibCore::library_namespace()+".LegendreGaussEndP4.Point");

ComponentBuilder<Point<LegendreGaussEnd,5>,mesh::ShapeFunction,LibCore>
PointP5_builder(LibCore::library_namespace()+".LegendreGaussEndP5.Point");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Line<LegendreGaussEnd,0>,mesh::ShapeFunction,LibCore>
LineP0_builder(LibCore::library_namespace()+".LegendreGaussEndP0.Line");

ComponentBuilder<Line<LegendreGaussEnd,1>,mesh::ShapeFunction,LibCore>
LineP1_builder(LibCore::library_namespace()+".LegendreGaussEndP1.Line");

ComponentBuilder<Line<LegendreGaussEnd,2>,mesh::ShapeFunction,LibCore>
LineP2_builder(LibCore::library_namespace()+".LegendreGaussEndP2.Line");

ComponentBuilder<Line<LegendreGaussEnd,3>,mesh::ShapeFunction,LibCore>
LineP3_builder(LibCore::library_namespace()+".LegendreGaussEndP3.Line");

ComponentBuilder<Line<LegendreGaussEnd,4>,mesh::ShapeFunction,LibCore>
LineP4_builder(LibCore::library_namespace()+".LegendreGaussEndP4.Line");

ComponentBuilder<Line<LegendreGaussEnd,5>,mesh::ShapeFunction,LibCore>
LineP5_builder(LibCore::library_namespace()+".LegendreGaussEndP5.Line");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Quad<LegendreGaussEnd,0>,mesh::ShapeFunction,LibCore>
QuadP0_builder(LibCore::library_namespace()+".LegendreGaussEndP0.Quad");

ComponentBuilder<Quad<LegendreGaussEnd,1>,mesh::ShapeFunction,LibCore>
QuadP1_builder(LibCore::library_namespace()+".LegendreGaussEndP1.Quad");

ComponentBuilder<Quad<LegendreGaussEnd,2>,mesh::ShapeFunction,LibCore>
QuadP2_builder(LibCore::library_namespace()+".LegendreGaussEndP2.Quad");

ComponentBuilder<Quad<LegendreGaussEnd,3>,mesh::ShapeFunction,LibCore>
QuadP3_builder(LibCore::library_namespace()+".LegendreGaussEndP3.Quad");

ComponentBuilder<Quad<LegendreGaussEnd,4>,mesh::ShapeFunction,LibCore>
QuadP4_builder(LibCore::library_namespace()+".LegendreGaussEndP4.Quad");

ComponentBuilder<Quad<LegendreGaussEnd,5>,mesh::ShapeFunction,LibCore>
QuadP5_builder(LibCore::library_namespace()+".LegendreGaussEndP5.Quad");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Hexa<LegendreGaussEnd,0>,mesh::ShapeFunction,LibCore>
HexaP0_builder(LibCore::library_namespace()+".LegendreGaussEndP0.Hexa");

ComponentBuilder<Hexa<LegendreGaussEnd,1>,mesh::ShapeFunction,LibCore>
HexaP1_builder(LibCore::library_namespace()+".LegendreGaussEndP1.Hexa");

ComponentBuilder<Hexa<LegendreGaussEnd,2>,mesh::ShapeFunction,LibCore>
HexaP2_builder(LibCore::library_namespace()+".LegendreGaussEndP2.Hexa");

ComponentBuilder<Hexa<LegendreGaussEnd,3>,mesh::ShapeFunction,LibCore>
HexaP3_builder(LibCore::library_namespace()+".LegendreGaussEndP3.Hexa");

ComponentBuilder<Hexa<LegendreGaussEnd,4>,mesh::ShapeFunction,LibCore>
HexaP4_builder(LibCore::library_namespace()+".LegendreGaussEndP4.Hexa");

ComponentBuilder<Hexa<LegendreGaussEnd,5>,mesh::ShapeFunction,LibCore>
HexaP5_builder(LibCore::library_namespace()+".LegendreGaussEndP5.Hexa");

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

