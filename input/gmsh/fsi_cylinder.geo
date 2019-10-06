// This is my first trial of arterial wall geometry
msize = 0.1;

Point(1) = {0, 0, 0, msize};
Point(2) = {1, 0, 0, msize};
Point(3) = {0, 1, 0, msize};
Point(4) = {-1, 0, 0, msize};
Point(5) = {0, -1, 0, msize};
Point(6) = {1.2, 0, 0, msize};
Point(7) = {0, 1.2, 0, msize};
Point(8) = {-1.2, 0, 0, msize};
Point(9) = {0, -1.2, 0, msize};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Line Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};
Line Loop(2) = {6, 7, 8, 5};
Plane Surface(2) = {1, 2};

Extrude {0, 0, 10} {
  Point{1}; Surface{1}; Surface{2}; 
}


Physical Surface("ftop") = {31};
Physical Surface("fbot") = {1};
Physical Surface("fwall") = {22, 18, 26, 30};
Physical Surface("stop") = {73};
Physical Surface("sbot") = {2};
Physical Surface("swall") = {64, 68, 72, 60};

Physical Volume("fluid") = {1};
Physical Volume("solid") = {2};
