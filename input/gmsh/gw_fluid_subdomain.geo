Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {-1, 0, 0};
Point(5) = {0, -1, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

h_bl = 0.8;

Point(10) = {h_bl, 0, 0, 1.0};
Point(11) = {0, h_bl, 0, 1.0};
Point(12) = {-h_bl, 0, 0, 1.0};
Point(13) = {0, -h_bl, 0, 1.0};

Circle(13) = {10, 1, 11};
Circle(14) = {11, 1, 12};
Circle(15) = {12, 1, 13};
Circle(16) = {13, 1, 10};

Line(17) = {10, 2};
Line(18) = {11, 3};
Line(19) = {12, 4};
Line(20) = {13, 5};

Transfinite Line {-17,-18,-19,-20} = 8 Using Progression 1.05;
Transfinite Line {1,2,3,4,5,6,7,8,13,14,15,16} = 31 Using Progression 1;


Line Loop(5) = {1, -18, -13, 17};
Plane Surface(5) = {5};
Line Loop(6) = {2, -19, -14, 18};
Plane Surface(6) = {6};
Line Loop(7) = {3, -20, -15, 19};
Plane Surface(7) = {7};
Line Loop(8) = {4, -17, -16, 20};
Plane Surface(8) = {8};
Line Loop(9) = {16, 13, 14, 15};
Plane Surface(9) = {9};

Transfinite Surface {5,6,7,8};

Extrude {0, 0, 10} {
  Surface{5,6,7,8,9}; Layers{150};
}

Physical Surface("wall") = {73, 51, 29, 95};
Physical Surface("outlet") = {130, 108, 42, 64, 86};
Physical Surface("inlet") = {7, 6, 9, 8, 5};
Physical Volume("vol") = {4, 5, 2, 1, 3};
