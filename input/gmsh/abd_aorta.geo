Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFactor = 0.1;

Merge "comp.stl";
CreateTopology;

//Extrude outward and inward with 4 layers of thickness 0.5
out1[] = Extrude{Surface{1}; Layers{2, 0.1}; Using Index[0]; };

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

Line Loop(400)={b1[2]};
Plane Surface(401)={400};

Line Loop(500)={b1[3]};
Plane Surface(501)={500};

Line Loop(600)={b1[4]};
Plane Surface(601)={600};

Line Loop(700)={b1[5]};
Plane Surface(701)={700};

Line Loop(750)={b1[6]};
Plane Surface(751)={750};

Line Loop(800)={b1[7]};
Plane Surface(801)={800};

Line Loop(850)={b1[8]};
Plane Surface(851)={850};


//create inside volume
Surface Loop(999)={1, 201,301,401, 501, 601, 701, 751, 801, 851};
Volume(1000) = {999};

//save only physicals
Physical Surface("inlet")={ 201 };
Physical Surface("outlet_1")={ 301 };
Physical Surface("outlet_2")={ 401 };
Physical Surface("outlet_3")={ 501 };
Physical Surface("outlet_4")={ 601 };
Physical Surface("outlet_5")={ 701 };
Physical Surface("outlet_6")={ 751 };
Physical Surface("outlet_7")={ 801 };
Physical Surface("outlet_8")={ 851 };
Physical Surface("fwall")={ 1 };
Physical Surface("swall")={ out1[0] };
Physical Surface("sbl")={ out1[2], out1[3], out1[4], out1[5], out1[6], out1[7], out1[8], out1[9], out1[10] }; 
Physical Volume("fluid")={1000};
Physical Volume("solid")={out1[1]};

//eof
