Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFactor = 0.05;

Merge "demo.stl";
CreateTopology;

//Extrude outward and inward with 4 layers of thickness 0.5
out1[] = Extrude{Surface{1}; Layers{2, 0.2}; Using Index[0]; };

b1[] = Boundary{ Surface{1}; };


For ii In { 0 : (#b1[] - 1) }
Printf("b1 number %g = %g", ii, b1[ii]);
EndFor


//create inlet faces
Line Loop(200)={b1[0]};
Plane Surface(201)={200};

Line Loop(300)={b1[1]};
Plane Surface(301)={300};

Line Loop(400)={b1[2]};
Plane Surface(401)={400};


//create inside volume
Surface Loop(999)={1, 201,301,401};
Volume(1000) = {999};

//save only physicals
Physical Surface("inlet")={ 201 };
Physical Surface("outlet")={ 301, 401 };
Physical Surface("sbl")={ out1[2], out1[3], out1[4] }; 
Physical Volume("fluid")={1000};
Physical Volume("solid")={out1[1]};

//eof
