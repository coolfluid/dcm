// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


/// @file LegendreGaussLobatto.cpp
///
/// This file creates and registers shape functions for tensorial
/// elements with Legendre-Gauss-Lobatto flux-points
///
/// Shape-functions are named e.g. "cf3.sdm.core.LegendreGaussLobattoP2"

#include "cf3/sdm/core/Tensorial.hpp"
#include "cf3/sdm/core/LegendreGaussLobatto.hpp"
#include "cf3/common/Builder.hpp"

using namespace cf3::common;

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Point<LegendreGaussLobatto,0>,mesh::ShapeFunction,LibCore>
PointP0_builder(LibCore::library_namespace()+".LegendreGaussLobattoP0.Point");

ComponentBuilder<Point<LegendreGaussLobatto,1>,mesh::ShapeFunction,LibCore>
PointP1_builder(LibCore::library_namespace()+".LegendreGaussLobattoP1.Point");

ComponentBuilder<Point<LegendreGaussLobatto,2>,mesh::ShapeFunction,LibCore>
PointP2_builder(LibCore::library_namespace()+".LegendreGaussLobattoP2.Point");

ComponentBuilder<Point<LegendreGaussLobatto,3>,mesh::ShapeFunction,LibCore>
PointP3_builder(LibCore::library_namespace()+".LegendreGaussLobattoP3.Point");

ComponentBuilder<Point<LegendreGaussLobatto,4>,mesh::ShapeFunction,LibCore>
PointP4_builder(LibCore::library_namespace()+".LegendreGaussLobattoP4.Point");

ComponentBuilder<Point<LegendreGaussLobatto,5>,mesh::ShapeFunction,LibCore>
PointP5_builder(LibCore::library_namespace()+".LegendreGaussLobattoP5.Point");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Line<LegendreGaussLobatto,0>,mesh::ShapeFunction,LibCore>
LineP0_builder(LibCore::library_namespace()+".LegendreGaussLobattoP0.Line");

ComponentBuilder<Line<LegendreGaussLobatto,1>,mesh::ShapeFunction,LibCore>
LineP1_builder(LibCore::library_namespace()+".LegendreGaussLobattoP1.Line");

ComponentBuilder<Line<LegendreGaussLobatto,2>,mesh::ShapeFunction,LibCore>
LineP2_builder(LibCore::library_namespace()+".LegendreGaussLobattoP2.Line");

ComponentBuilder<Line<LegendreGaussLobatto,3>,mesh::ShapeFunction,LibCore>
LineP3_builder(LibCore::library_namespace()+".LegendreGaussLobattoP3.Line");

ComponentBuilder<Line<LegendreGaussLobatto,4>,mesh::ShapeFunction,LibCore>
LineP4_builder(LibCore::library_namespace()+".LegendreGaussLobattoP4.Line");

ComponentBuilder<Line<LegendreGaussLobatto,5>,mesh::ShapeFunction,LibCore>
LineP5_builder(LibCore::library_namespace()+".LegendreGaussLobattoP5.Line");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Quad<LegendreGaussLobatto,0>,mesh::ShapeFunction,LibCore>
QuadP0_builder(LibCore::library_namespace()+".LegendreGaussLobattoP0.Quad");

ComponentBuilder<Quad<LegendreGaussLobatto,1>,mesh::ShapeFunction,LibCore>
QuadP1_builder(LibCore::library_namespace()+".LegendreGaussLobattoP1.Quad");

ComponentBuilder<Quad<LegendreGaussLobatto,2>,mesh::ShapeFunction,LibCore>
QuadP2_builder(LibCore::library_namespace()+".LegendreGaussLobattoP2.Quad");

ComponentBuilder<Quad<LegendreGaussLobatto,3>,mesh::ShapeFunction,LibCore>
QuadP3_builder(LibCore::library_namespace()+".LegendreGaussLobattoP3.Quad");

ComponentBuilder<Quad<LegendreGaussLobatto,4>,mesh::ShapeFunction,LibCore>
QuadP4_builder(LibCore::library_namespace()+".LegendreGaussLobattoP4.Quad");

ComponentBuilder<Quad<LegendreGaussLobatto,5>,mesh::ShapeFunction,LibCore>
QuadP5_builder(LibCore::library_namespace()+".LegendreGaussLobattoP5.Quad");

////////////////////////////////////////////////////////////////////////////////

ComponentBuilder<Hexa<LegendreGaussLobatto,0>,mesh::ShapeFunction,LibCore>
HexaP0_builder(LibCore::library_namespace()+".LegendreGaussLobattoP0.Hexa");

ComponentBuilder<Hexa<LegendreGaussLobatto,1>,mesh::ShapeFunction,LibCore>
HexaP1_builder(LibCore::library_namespace()+".LegendreGaussLobattoP1.Hexa");

ComponentBuilder<Hexa<LegendreGaussLobatto,2>,mesh::ShapeFunction,LibCore>
HexaP2_builder(LibCore::library_namespace()+".LegendreGaussLobattoP2.Hexa");

ComponentBuilder<Hexa<LegendreGaussLobatto,3>,mesh::ShapeFunction,LibCore>
HexaP3_builder(LibCore::library_namespace()+".LegendreGaussLobattoP3.Hexa");

ComponentBuilder<Hexa<LegendreGaussLobatto,4>,mesh::ShapeFunction,LibCore>
HexaP4_builder(LibCore::library_namespace()+".LegendreGaussLobattoP4.Hexa");

ComponentBuilder<Hexa<LegendreGaussLobatto,5>,mesh::ShapeFunction,LibCore>
HexaP5_builder(LibCore::library_namespace()+".LegendreGaussLobattoP5.Hexa");

////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

