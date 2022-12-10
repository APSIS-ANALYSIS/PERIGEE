// This is the gmsh geometry file for studying the turbulent flow
// passing a square cylinder in 3D.
// Reference: Turbulent flow around a square cylinder at Reynolds
//            number 22000: A DNS study. Computers and Fluids, 123 
//            (2015):87-98.
// Definition of mesh parameters
ee_WALL    = 100; // number of elements along the cylinder wall surface
ee_BL      = 64;  // number of boundary layers
ee_updown  = 91;  // number of elements in streamwise direction
ee_lateral = 101; // number of elements in cross-stream direction
ee_span    = 80;  // number of elements in spanwise direction

Point(1) = {-10.5, -27, -2};
Point(2) = {20.0, -27, -2};
Point(3) = {20.0, 27, -2};
Point(4) = {-10.5, 27, -2};

Point(5) = {-0.5, -0.5, -2};
Point(6) = {0.5, -0.5, -2};
Point(7) = {0.5, 0.5, -2};
Point(8) = {-0.5, 0.5, -2};

Line(5) = {7, 8};
Line(6) = {8, 5};
Line(7) = {5, 6};
Line(8) = {6, 7};

Line Loop(2) = {5, 6, 7, 8};

Point(9) = {-0.75, -0.75, -2};
Point(10) = {0.75, -0.75, -2};
Point(11) = {0.75, 0.75, -2};
Point(12) = {-0.75, 0.75, -2};

Line(9) = {11, 12};
Line(10) = {12, 9};
Line(11) = {9, 10};
Line(12) = {10, 11};

Line(13) = {4, 3};
Line(14) = {3, 2};
Line(15) = {2, 1};
Line(16) = {1, 4};
Line(17) = {5, 9};
Line(18) = {6, 10};
Line(19) = {7, 11};
Line(20) = {8, 12};

Transfinite Curve {7, 6, 5, 8} = ee_WALL + 1 Using Progression 1;
Transfinite Curve {11, 12, 10, 9, 29, 26, 35, 32} = ee_WALL+1 Using Progression 1;
Transfinite Curve {17, 18, 19, 20} = ee_BL+1 Using Progression 1.05;

Transfinite Curve {15, 13} = ee_updown Using Progression 1;

Transfinite Curve {14, 16} = ee_lateral Using Progression 1;

Curve Loop(3) = {17, 11, -18, -7};
Plane Surface(1) = {3};
Curve Loop(4) = {8, 19, -12, -18};
Plane Surface(2) = {4};
Curve Loop(5) = {5, 20, -9, -19};
Plane Surface(3) = {5};
Curve Loop(6) = {6, 17, -10, -20};
Plane Surface(4) = {6};
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Curve Loop(7) = {15, 16, 13, 14};
Curve Loop(8) = {11, 12, 9, 10};
Plane Surface(5) = {7, 8};

Extrude {0,0,4}{
  Surface{1, 2, 3, 4, 5}; Layers{ee_span}; 
}

Physical Surface("wall",   1) = {73, 51, 95, 41};
Physical Surface("inlet",  2) = {125};
Physical Surface("outlet", 3) = {133};
Physical Surface("top",    4) = {129};
Physical Surface("bottom", 5) = {121};
Physical Surface("front",  6) = {150, 108, 86, 64, 42};
Physical Surface("back",   7) = {5, 2, 3, 4, 1};
Physical Volume("fluid",   8) = {1, 2, 3, 4, 5};

// EOF
