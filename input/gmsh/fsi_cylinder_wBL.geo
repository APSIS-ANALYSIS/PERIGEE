Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {-1, 0, 0};
Point(5) = {0, -1, 0};
Point(6) = {1.2, 0, 0};
Point(7) = {0, 1.2, 0};
Point(8) = {-1.2, 0, 0};
Point(9) = {0, -1.2, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Line(9) = {7, 3};
Line(10) = {6, 2};
Line(11) = {9, 5};
Line(12) = {8, 4};

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

Transfinite Line {9,10,11,12} = 4 Using Progression 1;
Transfinite Line {-17,-18,-19,-20} = 6 Using Progression 1.05;
Transfinite Line {1,2,3,4,5,6,7,8,13,14,15,16} = 22 Using Progression 1;

Line Loop(1) = {5, 9, -1, -10};
Plane Surface(1) = {1};
Line Loop(2) = {9, 2, -12, -6};
Plane Surface(2) = {2};
Line Loop(3) = {12, 3, -11, -7};
Plane Surface(3) = {3};
Line Loop(4) = {11, 4, -10, -8};
Plane Surface(4) = {4};
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

Transfinite Surface {1,2,3,4,5,6,7,8};

Extrude {0, 0, 10} {
  Surface{1,2,3,4,5,6,7,8,9}; Layers{150};
}


Physical Surface("ftop") = {218, 152, 130, 196, 174};
Physical Surface("fbot") = {9, 5, 6, 7, 8};
Physical Surface("fwall") = {37, 99, 77, 55};
Physical Surface("stop") = {64, 86, 108, 42};
Physical Surface("sbot") = {2, 1, 4, 3};
Physical Surface("swall") = {29, 107, 85, 63};

Physical Volume("fluid") = {9, 5, 8, 7, 6};
Physical Volume("solid") = {1, 2, 3, 4};

// EOF
