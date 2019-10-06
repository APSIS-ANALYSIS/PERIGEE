// Gmsh for cylinder in a rectangular domain

Nx1 = 21; Rx1 = 1.00;
Nx2 = 41; Rx2 = 1.00;
Ny  = 21; Ry  = 5.00;
Nb  = 21; Rb  = 0.90;
Nc  = 21; Rc  = 1.00;


// Point for the external boundaries
Point(1) = {-10, -10, 0, 1.0};
Point(2) = {10, -10, 0, 1.0};
Point(3) = {40, -10, 0, 1.0};
Point(4) = {-10, 10, 0, 1.0};
Point(5) = {10, 10, 0, 1.0};
Point(6) = {40, 10, 0, 1.0};

// Point for the cylinder
cypt = 0.35355339;
Point(7) = {-cypt, -cypt, 0, 1.0};
Point(8) = {cypt, -cypt, 0, 1.0};
Point(9) = {-cypt, cypt, 0, 1.0};
Point(10) = {cypt, cypt, 0, 1.0};
Point(11) = {0, 0, 0, 1.0};

// Outer lines
Line(1) = {1, 2}; Transfinite Line {1} = Ny Using Bump Ry;
Line(2) = {2, 3}; Transfinite Line {2} = Nx2 Using Progression Rx2;
Line(3) = {4, 5}; Transfinite Line {3} = Ny Using Bump Ry;
Line(4) = {5, 6}; Transfinite Line {4} = Nx2 Using Progression Rx2;
Line(5) = {1, 4}; Transfinite Line {5} = Ny Using Bump Ry;
Line(6) = {2, 5}; Transfinite Line {6} = Ny Using Bump Ry;
Line(7) = {3, 6}; Transfinite Line {7} = Ny Using Bump Ry;

//Cylinder
Circle(8) = {8, 11, 10}; Transfinite Line {8} = Nc Using Progression Rc;
Circle(9) = {10, 11, 9}; Transfinite Line {9} = Nc Using Progression Rc;
Circle(10) = {9, 11, 7}; Transfinite Line {10} = Nc Using Progression Rc;
Circle(11) = {7, 11, 8}; Transfinite Line {11} = Nc Using Progression Rc;

//Block lines
Line(12) = {1, 7}; Transfinite Line {12} = Nb Using Progression Rb;
Line(13) = {2, 8}; Transfinite Line {13} = Nb Using Progression Rb;
Line(14) = {5, 10}; Transfinite Line {14} = Nb Using Progression Rb;
Line(15) = {4, 9}; Transfinite Line {15} = Nb Using Progression Rb;

// Define surfaces
Line Loop(1) = {12, 11, -13, -1};
Plane Surface(1) = {1};
Line Loop(2) = {6, 14, -8, -13};
Plane Surface(2) = {2};
Line Loop(3) = {14, 9, -15, 3};
Plane Surface(3) = {3};
Line Loop(4) = {15, 10, -12, 5};
Plane Surface(4) = {4};
Line Loop(5) = {6, 4, -7, -2};
Plane Surface(5) = {5};

Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {5};
//Recombine Surface {3, 4, 1, 2, 5};


Extrude {0, 0, 1} {
  Surface{3}; Surface{4}; Surface{1}; Surface{2}; Surface{5};
  Layers{2};
}

Physical Surface("inlet") = {58};
Physical Surface("outlet") = {120};
Physical Surface("solidwall") = {50, 98, 28, 72, 80, 124, 5, 125, 1, 81, 4, 59, 3, 37, 2, 103, 36, 116};
Physical Volume("vol1") = {1, 2, 3, 4};
Physical Volume("vol2") = {5};
