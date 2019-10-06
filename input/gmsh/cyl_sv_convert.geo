Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFactor = 0.1;

Merge "cyl.stl";
CreateTopology;

//Extrude outward and inward with 4 layers of thickness 0.5
out1[] = Extrude{Surface{1}; Layers{2, 0.2}; Using Index[0]; };

b1[] = Boundary{ Surface{1}; };

// Print stl info
For ii In { 0 : (#b1[] - 1) }
Printf("b1 number %g = %g", ii, b1[ii]);
EndFor

For ii In { 0 : (#out1[] - 1) }
Printf("out1 number %g = %g", ii, out1[ii]);
EndFor

//create inlet faces
Line Loop(200)={b1[0]};
Plane Surface(201)={200};

Line Loop(300)={b1[1]};
Plane Surface(301)={300};

//create inside volume
Surface Loop(999)={1,201,301};
Volume(1000) = {999};

//save only physicals
Physical Surface("inlet")={ 201 };
Physical Surface("outlet_1")={ 301 };
Physical Surface("fwall")={ 1 };
Physical Surface("swall")={ out1[0] };
Physical Surface("sbl")={ out1[2], out1[3] }; 
Physical Volume("fluid")={1000};
Physical Volume("solid")={out1[1]};

//eof
