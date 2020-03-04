Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.ElementOrder=1;
Mesh.CharacteristicLengthFactor = 0.5;

Point(1) = {0, 0, 0};
Point(2) = {2, 0, 0};
Point(3) = {-2, 0, 0};
Point(4) = {0, 2, 0};
Point(5) = {0, -2, 0};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};

Point(6) = {1.8, 0, 0};
Point(7) = {0, 1.8, 0};
Point(8) = {-1.8, 0, 0};
Point(9) = {0, -1.8, 0};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Line(9)  = {4, 7};
Line(10) = {2, 6};
Line(11) = {5, 9};
Line(12) = {3, 8};

// Make 4 radial lines transfinite
Transfinite Line {9,10,11,12} = 5 Using Progression 1.25;

Line Loop(1) = {5, 6, 7, 8};

Plane Surface(1) = {1};

Line Loop(2) = {1, 9, -5, -10};

Plane Surface(2) = {2};

Line Loop(3) = {2, 12, -6, -9};

Plane Surface(3) = {3};

Line Loop(4) = {-12, 3, 11, -7};

Plane Surface(4) = {4};

Line Loop(5) = {-11, 4, 10, -8};

Plane Surface(5) = {5};

// Make the annular region transfinite so we get structured mesh here
Transfinite Surface {2,3,4,5};

// Make the two circle transfinite
Transfinite Line {1,2,3,4,5,6,7,8} = 10 Using Progression 1;

Extrude {0, 0, 30} {
  Surface{1,2,3,4,5}; Layers{150};
}

Physical Surface("wall") = {65, 43, 113, 91};
Physical Surface("inlet") = {34, 78, 56, 122, 100};
Physical Surface("outlet") = {1, 5, 2, 3, 4};
Physical Volume("cyl") = {1, 4, 5, 2, 3};

// EOF
