Point(1) = {0,0,0};
Point(2) = {1, 0, 0};
Point(3) = {-1, 0, 0};
Point(4) = {0, 1, 0};
Point(5) = {0, -1, 0};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};


Line Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};


// boundary layer
Field[4] = BoundaryLayer;
Field[4].EdgesList = {1,2,3,4};  // FacesList  for 3D
Field[4].hwall_n = 0.02;
Field[4].hwall_t = 0.01;
Field[4].thickness = 0.1;
Field[4].ratio = 1.2;
 
BoundaryLayer Field = 4;

