gs = 0.1;
Point(1) = {0,0,0, gs};
Point(2) = {1,0,0, gs};
Point(3) = {1,1,0, gs};
Point(4) = {0,1,0, gs};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Point(5) = {2.0, 1.0, 0, gs};
Point(6) = {2.0, 0.0, 0, gs};
Point(7) = {2.0, 0.0, 0, gs};
Line(9) = {3, 5};
Line(10) = {5, 6};
Line(11) = {6, 2};

// Define two square domain
Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};
Line Loop(2) = {6, 9, 10, 11};
Plane Surface(2) = {2};

Transfinite Line {8, 5, 6, 7} = 16 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("face1") = {1};
Physical Surface("face2") = {2};
