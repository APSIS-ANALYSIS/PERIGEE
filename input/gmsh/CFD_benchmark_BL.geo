Mesh.Algorithm3D=4;
Mesh.Optimize=1;
Mesh.CharacteristicLengthFactor=0.05;
Mesh.ElementOrder=1;

//Mesh.Smoothing = 100;
//Mesh.OptimizeNetgen=1;

// Point mesh size at the center in the expanded tube
//pt_len_out = 0.01;

// Boundary layer progression value and number of layers
val_prog = 1.2; 
num_bl = 8;
bl_cirm_num = 5; // number of elements in BL circumferential direction
bl_thickness_ratio = 0.15;

out_val_prog=1.5; 
out_num_bl=6;
out_bl_cirm_num=10;
out_bl_thickness_ratio = 0.1;

// Axial extrusion layer lenth
len = 0.1;

// Geometry Parameter
delta = -6.2685;

df=0.6;
di=0.6;
dt=0.2;

// Boundary layer thickness in inflow and throat
bli  = 0.6 * bl_thickness_ratio;
blt  = 0.2 * bl_thickness_ratio;
obli = 0.6 * out_bl_thickness_ratio;

l1=20*di;
l2=2.2685;
l3=4;
l4=30*di;

num_layer_l1 = l1 / len;
num_layer_l2 = l2 / len;
num_layer_l3 = l3 / len;
num_layer_l4 = l4 / len;

Point(1) = {0,  0, 0 + delta};
Point(2) = {di-bli, 0, 0 + delta};
Point(6) = {0,  0, l2 + delta};
Point(7) = {dt-blt, 0, l2 + delta};
Line(1) = {1,2};
Line(2) = {1,6};
Line(3) = {6,7};
Line(4) = {7,2};
Line Loop(1) = {4, -1, 2, 3};
Plane Surface(1) = {1};

Point(8) = {di, 0, 0 + delta};
Point(9) = {dt, 0, l2 + delta};

Line(22) = {8, 2};
Line(23) = {9, 7};
Line(24) = {9, 8};
Line Loop(25) = {24, 22, -4, -23};
Plane Surface(26) = {25};

Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{1}; Surface{26}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{43}; Surface{65}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{82}; Surface{104}; }
Extrude {{0,0,1}, {0,0,0+delta}, Pi/2} { Surface{121}; Surface{143}; }

// Make BL mesh
Transfinite Line {22,124,85,46, -48,23,-87,-126} = num_bl Using Progression val_prog;
Transfinite Line {34,51,73,90,112,129,151,167} = bl_cirm_num Using Progression 1;
Transfinite Line {72,89,33,50,150,166,111,128} = bl_cirm_num Using Progression 1;
Transfinite Surface {56,95,134,172,64,180,142,103};

// Inflow region
Extrude {0, 0, -l1} {
  Surface{77}; Surface{95}; Surface{116}; Surface{134};
Surface{38}; Surface{56}; Surface{155}; Surface{172}; Layers{num_layer_l1};
}

// Throat
Extrude {0, 0, l3} {
  Surface{81}; Surface{42}; Surface{120}; Surface{159};
Surface{64}; Surface{103}; Surface{142}; Surface{180}; Layers{num_layer_l3};
}

// Sudden expansion
Point(98)  = {di-obli, 0, 0};
Point(99)  = {obli-di, 0, 0};
Point(100) = {0, di-obli, 0};
Point(101) = {0, obli-di, 0};

Circle(474) = {100, 61, 99};
Circle(475) = {99, 61, 101};
Circle(476) = {101, 61, 98};
Circle(477) = {98, 61, 100};

Point(102) = {di, 0, 0};
Point(103) = {-di, 0, 0};
Point(104) = {0, di, 0};
Point(105) = {0, -di, 0};

Circle(478) = {104, 61, 103};
Circle(479) = {103, 61, 105};
Circle(480) = {105, 61, 102};
Circle(481) = {102, 61, 104};

Line(482) = {104, 100};
Line(483) = {103, 99};
Line(484) = {105, 101};
Line(485) = {102, 98};

Line Loop(29) = {478, 483, -474, -482};
Plane Surface(495) = {29};
Line Loop(30) = {483, 475, -484, -479};
Plane Surface(496) = {30};
Line Loop(31) = {484, 476, -485, -480};
Plane Surface(497) = {31};
Line Loop(32) = {485, 477, -482, -481};
Plane Surface(498) = {32};

Line Loop(33) = {474, 475, 476, 477};
Line Loop(34) = {473, 407, 429, 451};
Plane Surface(499) = {33, 34};

// After expansion cylinder BL definition
Transfinite Line {482,483,484,485} = out_num_bl Using Progression out_val_prog;
Transfinite Line {474,475,476,477,478,479,480,481} = out_bl_cirm_num Using Progression 1;
Transfinite Surface {495,496,497,498};

// Expansion volume by extrusion
Extrude {0, 0, l4} {
  Surface{387,353,370,404,470,492,426,448,495,496,497,498,499}; Layers{num_layer_l4};
}

//Characteristic Length {61} = pt_len_out;
//Characteristic Length {106} = pt_len_out;

// Define physical domains
Physical Surface("wall") = {698,664,742,720,499,497,498,496,495,483,417,439,461,130,91,52,168,218,257,335,296};
Physical Surface("outflow") = {677,699,721,743,785,655,589,611,633,550,533,516,567};
Physical Surface("inflow") = {258,219,297,336,236,314,275,197};
Physical Volume("vol") = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};


// EOF
