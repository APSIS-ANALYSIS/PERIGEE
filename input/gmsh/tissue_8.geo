Mesh.Algorithm=6;
Mesh.Optimize=1;
Mesh.ElementOrder=1;

Mesh.Smoothing = 100;
Mesh.OptimizeNetgen=1;
//Mesh.OptimizeThreshold =0.2;

ex = 40;
ey = ex * 6;
ez = ex * 20;

Point(1) = {0,0,0};
Point(2) = {0.15,0,0};
Point(3) = {0.15,0.025,0};
Point(4) = {0,0.025,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};

Transfinite Line {5, 7} = ey+1 Using Progression 1;
Transfinite Line {8, 6} = ex+1 Using Progression 1;
Transfinite Surface {1};

Extrude {0, 0, 0.5} {Surface{1}; Layers{ez};}

Physical Surface("top") = {30};
Physical Surface("bot") = {1};
Physical Surface("lef") = {21};
Physical Surface("rig") = {29};
Physical Surface("fro") = {25};
Physical Surface("bac") = {17};

Physical Volume("vol") = {1};

// EOF
