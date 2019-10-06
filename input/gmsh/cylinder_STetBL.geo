Point(1) = {0,0,0};
Point(2) = {1, 0, 0};
Point(3) = {-1, 0, 0};
Point(4) = {0, 1, 0};
Point(5) = {0, -1, 0};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};


Point(6) = {0.8, 0, 0};
Point(7) = {0, 0.8, 0};
Point(8) = {-0.8, 0, 0};
Point(9) = {0, -0.8, 0};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};


Line(9) = {4, 7};

Line(10) = {2, 6};

Line(11) = {5, 9};

Line(12) = {3, 8};


// Make 4 radial lines transfinite
Transfinite Line {9,10,11,12} = 5 Using Progression 1.2;


Line Loop(1) = {5, 6, 7, 8};

Plane Surface(1) = {1};

Line Loop(2) = {1, 9, -5, -10};

Plane Surface(2) = {2};

Line Loop(3) = {2, 12, -6, -9};

Plane Surface(3) = {3};

Line Loop(4) = {-12, 3, 11, -7};

Plane Surface(4) = {4};

Line Loop(5) = {-11, 4, 10, -8};

Plane Surface(5) = {5};

// Make the annular region transfinite so we get structured mesh here
Transfinite Surface {2,3,4,5};

// Make the two circle transfinite
Transfinite Line {1,2,3,4,5,6,7,8} = 20 Using Progression 1;


Extrude {0, 0, 1} {
  Surface{1,2,3,4,5}; Layers{15};
}



//Field[1] = Attractor;
//Field[1].NodesList = {6,7,8,9};
//Field[1].NNodesByEdge = 500;
//Field[1].EdgesList = {5,6,7,8};
//Field[1].dMax = 0.1;
//Field[1].dMin = 0.01;

//Field[2] = Threshold;
//Field[2].IField = 1;
//Field[2].LcMin = lc / 20;
//Field[2].LcMax = lc;
//Field[2].DistMin = 0.15;
//Field[2].DistMax = 0.5;

//Compound Surface(3) = {1, 2};


/*
// boundary layer
Field[4] = BoundaryLayer;
Field[4].EdgesList = {1,2,3,4};  // FacesList  for 3D
Field[4].hwall_n = 0.01;
Field[4].hwall_t = 0.01;
Field[4].thickness = 0.15;
Field[4].ratio = 1.1;
Field[4].Quads = 0;

BoundaryLayer Field = 4;
*/
/*
Field[1] = Attractor;
Field[1].EdgesList = {1,2,3,4};
//Field[1].NNodesByEdge = 50;
//Field[1].dMax = 0.3;
//Field[1].dMin = 0.1;

lc = 0.1;
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 10;
Field[2].LcMax = lc;
Field[2].DistMin = 0.01;
Field[2].DistMax = 0.1;

Background Field = 2;
*/



