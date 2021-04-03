ms = 0.35;

Point(1) = {0, 0, 0, ms};
Point(2) = {2, 0, 0, ms};
Point(3) = {-2, 0, 0, ms};
Point(4) = {0, 2, 0, ms};
Point(5) = {0, -2, 0, ms};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};


//Line Loop(1) = {2, 3, 4, 1};
//Plane Surface(1) = {1};

//Extrude {0, 0, 30.0} {Surface{1};}

//Physical Surface("inlet") = {26};
//Physical Surface("wall") = {21, 25, 17, 13};
//Physical Surface("outlet") = {1};
//Physical Volume("cyl") = {1};
//+
Line(5) = {1, 3};
Line(6) = {1, 5};
Line(7) = {1, 2};
Line(8) = {1, 4};

Curve Loop(1) = {5, -2, -8};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 3, -6};
Plane Surface(2) = {2};
Curve Loop(3) = {6, 4, -7};
Plane Surface(3) = {3};
Curve Loop(4) = {7, 1, -8};
Plane Surface(4) = {4};
Transfinite Curve {5, 8, 6, 7} = 16 Using Progression 1.0;
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

Extrude {0, 0, 30.0} {Surface{1,2,3,4};}

Physical Surface("inlet") = {25, 76, 59, 42};
Physical Surface("wall") = {71, 20, 37, 54};
Physical Surface("outlet") += {1, 4, 3, 2};
Physical Volume("cyl") = {4, 1, 2, 3};
