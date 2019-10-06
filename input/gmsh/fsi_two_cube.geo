ms = 0.5;
Point(1) = {0, 0, 0, ms};
Point(2) = {1, 0, 0, ms};
Point(3) = {1, 1, 0, ms};
Point(4) = {0, 1, 0, ms};

Point(5) = {0, 0, 1, ms};
Point(6) = {1, 0, 1, ms};
Point(7) = {1, 1, 1, ms};
Point(8) = {0, 1, 1, ms};

Point(9) = {0, 0, 2, ms};
Point(10) = {1, 0, 2, ms};
Point(11) = {1, 1, 2, ms};
Point(12) = {0, 1, 2, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 5};
Line(8) = {5, 6};
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 9};
Line(12) = {9, 10};
Line(13) = {1, 5};
Line(14) = {4, 8};
Line(15) = {3, 7};
Line(16) = {2, 6};
Line(17) = {7, 11};
Line(18) = {6, 10};
Line(19) = {8, 12};
Line(20) = {5, 9};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};


Line Loop(2) = {1, 16, -8, -13};
Plane Surface(2) = {2};
Line Loop(3) = {2, 15, -5, -16};
Plane Surface(3) = {3};
Line Loop(4) = {3, 14, -6, -15};
Plane Surface(4) = {4};
Line Loop(5) = {4, 13, -7, -14};
Plane Surface(5) = {5};
Line Loop(6) = {6, 7, 8, 5};
Plane Surface(6) = {6};
Line Loop(7) = {19, 11, -20, -7};
Plane Surface(7) = {7};
Line Loop(8) = {19, -10, -17, 6};
Plane Surface(8) = {8};
Line Loop(9) = {17, -9, -18, 5};
Plane Surface(9) = {9};
Line Loop(10) = {18, -12, -20, 8};
Plane Surface(10) = {10};
Line Loop(11) = {9, 10, 11, 12};
Plane Surface(11) = {11};


Surface Loop(1) = {10, 9, 8, 7, 11, 6};
Volume(1) = {1};
Surface Loop(2) = {2, 1, 3, 4, 5, 6};
Volume(2) = {2};

Physical Surface("ftop") = {11};
Physical Surface("fbot") = {10, 9, 8, 7};
Physical Surface("fwall") = {6};

Physical Surface("sbot") = {1};
Physical Surface("stop") = {2};
Physical Surface("swall") = {3, 4, 5};

Physical Volume("fluid") = {1};
Physical Volume("solid") = {2};
