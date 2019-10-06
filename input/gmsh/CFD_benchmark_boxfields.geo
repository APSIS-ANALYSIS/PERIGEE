Mesh.Algorithm3D=1;
Mesh.Optimize=1;
//Mesh.CharacteristicLengthFactor=0.05;
Mesh.ElementOrder=1;

//Mesh.Smoothing = 100;
//Mesh.OptimizeNetgen=1;

//ms = 0.2; // Entrance exit point size
//mse = 0.12; // Sudden expansion exterior size
//mm = 0.08;  // Converging size
//mt = 0.06; // Throad size

delta = -6.2685;

df=0.6;
di=0.6;
dt=0.2;

l1=20*di;
l2=2.2685;
l3=4;
l4=30*di;

Point(1) = {0, 0,0 + delta};
Point(2) = {di,0,0 + delta};
Point(6) = {0, 0,l2 + delta};
Point(7) = {dt,0,l2 + delta};
Line(1) = {1,2};
Line(2) = {1,6};
Line(3) = {6,7};
Line(4) = {7,2};
Line Loop(1) = {4, -1, 2, 3};
Plane Surface(1) = {1};

Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{1}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{21}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{38}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{55}; }

// Sudden expansion
Point(101) = {df,0,l2+l3+delta};
Point(102) = {0,df,l2+l3+delta};
Point(103) = {-df,0,l2+l3+delta};
Point(104) = {0,-df,l2+l3+delta};

Line Loop(2) = {28, 45, 62, 11};

Plane Surface(72) = {2};

Extrude {0, 0, l3} {
  Surface{54}; Surface{71}; Surface{20}; Surface{37};
}

Extrude {0, 0, -l1} {
  Surface{16}; Surface{33}; Surface{50}; Surface{67};
}

Circle(195) = {104, 105, 101};
Circle(196) = {101, 105, 102};
Circle(197) = {102, 105, 103};
Circle(198) = {103, 105, 104};

Line Loop(3) = {195, 196, 197, 198};
Line Loop(4) = {109, 126, 75, 92};

Plane Surface(209) = {3, 4};

Extrude {0, 0, l4} {
  Surface{209}; Surface{106}; Surface{123}; Surface{140}; Surface{89};
}

//Characteristic Length {118,105,126,106,111,6,7,8,12,10} = mt;
//Characteristic Length {1,2,9,11,13} = mm;
//Characteristic Length {126,127,132,136,134} = ms;
//Characteristic Length {144,139,137,149,138} = ms;

Field[1] = Box;
Field[1].VIn = 0.06;
Field[1].VOut = 0.1;
Field[1].XMax = 0.7;
Field[1].XMin = -0.7;
Field[1].YMax = 0.7;
Field[1].YMin = -0.7;
Field[1].ZMax = 1;
Field[1].ZMin = -7;

Field[2] = Box;
Field[2].VIn = 0.03;
Field[2].VOut = 0.1;
Field[2].XMax = 0.7;
Field[2].XMin = -0.7;
Field[2].YMax = 0.7;
Field[2].YMin = -0.7;
Field[2].ZMax = 0;
Field[2].ZMin = -4;

Field[7] = Min;
Field[7].FieldsList = {1,2};
Background Field = 7;

// Define physical domains
Physical Surface("wall") = {152, 203, 186, 169, 30, 47, 64, 13, 118, 135, 84, 101, 209, 222, 226, 230, 234};
Physical Surface("outflow") = {251, 285, 302, 319, 268};
Physical Surface("inflow") = {157, 174, 191, 208};
Physical Volume("vol") = {10, 9, 12, 11, 2, 1, 7, 8, 3, 5, 6, 4, 15, 16, 17, 13, 14};


