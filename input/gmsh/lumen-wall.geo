Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.05;

Merge "wall.stl";
Merge "inlet.stl";
Merge "outlet.stl";

// Remove duplicate vertices (will make surface mesh watertight)
Coherence Mesh;

// Refine by subdivision
//RefineMesh;

// Create volume
Surface Loop(1) = {1,2,3};
Volume(1) = {1};

Physical Surface("fbot") = {2};
Physical Surface("ftop") = {3};
Physical Surface("fwall") = {1};
Physical Volume("fluid") = {1};

// EOF
