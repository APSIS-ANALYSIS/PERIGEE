Mesh.Algorithm=5;
Mesh.Optimize=1;
Mesh.CharacteristicLengthFactor=0.1;

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

// Define two square domain
Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};

Physical Surface("face") = {1};
Physical Line("lef")     = {8};
Physical Line("top")     = {7};
Physical Line("rig")     = {6};
Physical Line("bot")     = {5};
