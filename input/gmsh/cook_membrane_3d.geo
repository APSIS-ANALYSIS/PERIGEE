Mesh.Algorithm=5;
Mesh.Optimize=1;
Mesh.CharacteristicLengthFactor=1.0;

Point(1) = {0,0,0};
Point(2) = {0.48,0.44,0};
Point(3) = {0.48,0.60,0};
Point(4) = {0,0.44,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};

Extrude {0, 0, 0.02} {Surface{1}; Layers{1};}

Physical Surface("top") = {30};
Physical Surface("bot") = {1};
Physical Surface("lef") = {17};
Physical Surface("rig") = {25};
Physical Surface("fro") = {21};
Physical Surface("bac") = {29};

Physical Volume("vol") = {1};

// EOF
