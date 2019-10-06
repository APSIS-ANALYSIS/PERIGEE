Mesh.Algorithm=5;
Mesh.Optimize=1;
Mesh.CharacteristicLengthFactor=1.0;
Mesh.ElementOrder=2;

ee = 50;

Point(1) = {0,0,0};
Point(2) = {1.0,0,0};
Point(3) = {1.0,1.0,0};
Point(4) = {0,1.0,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};

Transfinite Line {8, 5, 6, 7} = ee+1 Using Progression 1;
Transfinite Surface {1};

Extrude {0, 0, 1.0} {Surface{1}; Layers{ee};}

Physical Surface("top") = {30};
Physical Surface("bot") = {1};
Physical Surface("lef") = {21};
Physical Surface("rig") = {29};
Physical Surface("fro") = {25};
Physical Surface("bac") = {17};

Physical Volume("vol") = {1};

// EOF
