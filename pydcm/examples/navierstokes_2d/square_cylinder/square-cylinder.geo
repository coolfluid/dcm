DefineConstant [ lSquare = {1.0, Name "Dimensions/Square"} ];
DefineConstant [ xMin = {8.0, Name "Dimensions/before"} ];
DefineConstant [ xMax = {20.0, Name "Dimensions/after"} ];
DefineConstant [ yMax = {8.0, Name "Dimensions/top-bottom"} ];
DefineConstant [
  showLines = {Geometry.Lines, Choices {0,1}, Name "Options/Show lines",
              GmshOption "Geometry.Lines", AutoCheck 0}
  showPoints = {Geometry.Points, Choices {0,1}, Name "Options/Show points",
               GmshOption "Geometry.Points", AutoCheck 0}
  showNodes = {Mesh.Points, Choices {0,1}, Name "Options/Show nodes",
                GmshOption "Mesh.Points", AutoCheck 0}

  recombineAll = { 0 , Choices {0,1}, Name "Meshing/Recombine",
                GmshOption "Mesh.RecombineAll", AutoCheck 1}
  lcSquare = {0.1, Name "Meshing/Square"}
  lcInt1 = {0.5, Name "Meshing/Interior1"}
  lcInt2 = {1., Name "Meshing/Interior2"}
  lcExt = {4., Name "Meshing/Exterior"}
  blayer_thickness = {0.25, Name "Meshing/Blayer/Thickness"}
  blayer_ratio = {1.25, Name "Meshing/Blayer/Ratio"}
  blayer_hwall_n = {0.005, Name "Meshing/Blayer/hwall_n"}
  smoothing = {1 , Name "Meshing/Smoothing"}
];


Mesh.CharacteristicLengthExtendFromBoundary = 1;

// Meshing options
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8; // (0) Automatic (1) MeshAdapt (5) Delauney (6) Frontal (8) Delauney Quad (9) Packing of parallelograms
Mesh.RemeshAlgorithm = 1; // (0) no split (1) automatic (2) automatic only with metis
Mesh.RemeshParametrization = 0; // (0) harmonic (1) conformal spectral (7) conformal finite element
Mesh.RecombinationAlgorithm = 0;  // (0) standard (1) blossom
Mesh.SubdivisionAlgorithm = 1; // (0) None (1) All Quads (2) All Hexas
Mesh.RecombineAll=recombineAll;
Mesh.Smoothing = smoothing;

// Grid boundary layer size = d2-d1
d1 = lSquare/2.;
d2 = d1/5.*7.;

// Use to define 4 corners of the domain
xMax = xMax;
xMin = -xMin;
yMax = yMax;

// Define points, lines, line loops and surface for the square cylinder
p1 = newp; Point(p1) = {-d1,d1,0.0,lcSquare};
p2 = newp; Point(p2) = {-d1,-d1,0.0,lcSquare};
p3 = newp; Point(p3) = {d1,-d1,0.0,lcSquare};
p4 = newp; Point(p4) = {d1,d1,0.0,lcSquare};


// Define points, lines, line loops and surfaces for the intermediate and far field
p5 = newp; Point(p5) = {xMin,-yMax,0.0,lcExt};
p6 = newp; Point(p6) = {xMax,-yMax,0.0,lcExt};
p7 = newp; Point(p7) = {xMax,yMax,0.0,lcExt};
p8 = newp; Point(p8) = {xMin,yMax,0.0,lcExt};


l1 = newl; Line(l1) = {p1,p2}; line_1_2 = l1; line_2_1 = -l1;
l2 = newl; Line(l2) = {p2,p3}; line_2_3 = l2; line_3_1 = -l2;
l3 = newl; Line(l3) = {p3,p4}; line_3_4 = l3; line_4_3 = -l3;
l4 = newl; Line(l4) = {p4,p1}; line_4_1 = l4; line_1_4 = -l4;

l5 = newl; Line(l5) = {p5,p6}; line_5_6 = l5; line_6_5 = -l5;
l6 = newl; Line(l6) = {p6,p7}; line_6_7 = l6; line_7_6 = -l6;
l7 = newl; Line(l7) = {p7,p8}; line_7_8 = l7; line_8_7 = -l7;
l8 = newl; Line(l8) = {p8,p5}; line_8_5 = l8; line_5_8 = -l8;

ll1 = newl; Line Loop(ll1) = {line_1_2,line_2_3,line_3_4,line_4_1};
ll2 = newl; Line Loop(ll2) = {line_5_6,line_6_7,line_7_8,line_8_5};

surf = news; Plane Surface(surf) = {ll1,-ll2};

Transfinite Line {line_1_2,line_2_3,line_3_4,line_4_1} = lSquare/lcSquare Using Bump 0.2;

Field[1] = BoundaryLayer;
Field[1].EdgesList = {line_1_2,line_2_3,line_3_4,line_4_1};
Field[1].NodesList = {1,2,3,4};
Field[1].hfar = lcExt;
Field[1].hwall_n = blayer_hwall_n;
Field[1].hwall_t = blayer_hwall_n;
Field[1].thickness = blayer_thickness;
Field[1].ratio = blayer_ratio;
Field[1].fan_angle = 60;
Field[1].Quads = 1;
BoundaryLayer Field = 1;

// Define physical regions
Physical Line("farfield") = {line_5_6,line_7_8,line_8_5}; // Far field
Physical Line("outlet") = {line_6_7}; // Outlet
Physical Line("cylinder")  = {line_1_2, line_2_3, line_3_4, line_4_1}; // Cylinder
Physical Surface("domain") = {surf};

Field[2] = Box;
Field[2].VIn  = lcInt1;
Field[2].VOut = lcExt;
Field[2].XMax = 5*lSquare;
Field[2].XMin = -lSquare;
Field[2].YMax = 1.5*lSquare;
Field[2].YMin = -1.5*lSquare;

Field[3] = Box;
Field[3].VIn  = lcInt2;
Field[3].VOut = lcExt;
Field[3].XMax = 10*lSquare;
Field[3].XMin = -2.5*lSquare;
Field[3].YMax =  4*lSquare;
Field[3].YMin = -4*lSquare;

Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;
