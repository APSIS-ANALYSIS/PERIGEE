ms = 0.35;

Point(1) = {0, 0, 0, ms};
Point(2) = {Sqrt(2), 0, -Sqrt(2), ms};
Point(3) = {-Sqrt(2), 0, Sqrt(2), ms};
Point(4) = {Sqrt(6)/2, 1, Sqrt(6)/2, ms};
Point(5) = {-Sqrt(6)/2, -1, -Sqrt(6)/2, ms};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};

Extrude {Sqrt(2), -2*Sqrt(3), Sqrt(2)} {Surface{1};}

Physical Surface("inlet") = {26};
Physical Surface("wall") = {21, 25, 17, 13};
Physical Surface("outlet") = {1};
Physical Volume("cyl") = {1};


