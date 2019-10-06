//Mesh.Algorithm=6;
//Mesh.CharacteristicLengthMax=6.0e-2;
//Mesh.CharacteristicLengthMin=6.0e-2;

Point(1) = {-5, -6, 0};
Point(2) = {14.5, -6, 0};
Point(3) = {14.5, 6, 0};
Point(4) = {-5, 6, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Point(5) = {-0.5, -0.5, 0};
Point(6) = {0.5, -0.5, 0};
Point(7) = {0.5, 0.5, 0};
Point(8) = {-0.5, 0.5, 0};

Point(9) = {0.5, 0.03, 0};
Point(10) = {0.5, -0.03, 0};
Point(11) = {4.5, 0.03, 0};
Point(12) = {4.5, -0.03, 0};

Line(5) = {7, 8};
Line(6) = {8, 5};
Line(7) = {5, 6};
Line(8) = {6, 10};
Line(9) = {10, 9};
Line(10) = {9, 7};
Line(11) = {9, 11};
Line(12) = {11, 12};
Line(13) = {12, 10};
Line Loop(1) = {4, 1, 2, 3};
Line Loop(2) = {6, 7, 8, -13, -12, -11, 10, 5};
Plane Surface(1) = {1, 2};

Line Loop(3) = {13, 9, 11, 12};
Plane Surface(2) = {3};

Extrude {0,0,1}{
  Surface{1,2}; 
}

// Field 1 set the fluid mesh near the beam
Field[1] = Box;
Field[1].XMin = -0.8;
Field[1].XMax = 6;
Field[1].YMin = -0.8;
Field[1].YMax = 0.8;
Field[1].ZMin = 0;
Field[1].ZMax = 1;
Field[1].VIn = 0.1;
Field[1].VOut = 0.3;

Field[2] = Box;
Field[2].XMin = 0.5;
Field[2].XMax = 4.5;
Field[2].YMin = -0.03;
Field[2].YMax = 0.03;
Field[2].ZMin = 0;
Field[2].ZMax = 1;
Field[2].VIn = 0.06;
Field[2].VOut = 0.3;

Field[5] = Min;
Field[5].FieldsList = {1,2};

Background Field = 5;

Physical Surface("ftop") = {42};
Physical Surface("fbot") = {34};
Physical Surface("ffro") = {75};
Physical Surface("fbac") = {1};
Physical Surface("flef") = {30};
Physical Surface("frig") = {38};
Physical Surface("cube_wall") = {46, 74, 70, 66, 50};
Physical Surface("slef") = {88};
Physical Surface("srig") = {58};
Physical Surface("sfro") = {97};
Physical Surface("sbac") = {2};
Physical Surface("stop") = {54};
Physical Surface("sbot") = {62};

Surface Loop(1) = {75, 30, 1, 34, 38, 42, 46, 50, 74, 70, 66, 54, 58, 62};
Volume(3) = {1};
Surface Loop(2) = {62, 58, 54, 97, 88, 2};
Volume(4) = {2};


Physical Volume("fluid") = {1};
Physical Volume("solid") = {2};

// EOF
