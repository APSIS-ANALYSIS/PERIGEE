// This is a debug routine that generate four vol domains
msize = 0.1;

Point(1) = {0, 0, 0, msize};
Point(2) = {1, 0, 0, msize};
Point(3) = {1, 1, 0, msize};
Point(4) = {0, 1, 0, msize};

Point(5) = {2, 0, 0, msize};
Point(6) = {2, 1, 0, msize};

Point(7) = {0, 2, 0, msize};
Point(8) = {1, 2, 0, msize};


Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {2, 3};
Line(5) = {1, 4};
Line(6) = {4, 7};
Line(7) = {3, 8};
Line(8) = {4, 3};
Line(9) = {3, 6};
Line(10) = {7, 8};
Line(11) = {8, 6};

Line Loop(1) = {6, 10, -7, -8};
Plane Surface(1) = {1};

Line Loop(2) = {7, 11, -9};
Plane Surface(2) = {2};

Line Loop(3) = {5, 8, -4, -1};
Plane Surface(3) = {3};

Line Loop(4) = {4, 9, -3, -2};

Plane Surface(4) = {4};

Extrude {0, 0, 1} {
  Surface{1}; Surface{3}; Surface{4}; Surface{2}; 
}


Physical Surface("top") = {33, 94, 77, 55};

Physical Surface("bot") = {1, 2, 4, 3};

Physical Surface("lat") = {24, 89, 72, 76, 54, 42, 20};

Physical Volume("a") = {1};
Physical Volume("b") = {2};
Physical Volume("c") = {4};
Physical Volume("d") = {3};
