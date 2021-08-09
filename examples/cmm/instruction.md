# Instruction
This folder contains a suite of codes for performing patient-specific vascular FSI simulations using a small-strain membrane wall model.

## Pre-processing
In the preprocessor (source code is preprocess_sv_tets.cpp), one needs to prepare the mesh, assign proper boundary conditions, and do mesh partitioning. There are a few critical input arugments.
* `-cpu_size` determines the number of CPUs you want to run for the analysis solver.
* `-elem_type` determines the type of mesh you want to use.
* `-num_outlet` determines the number of outlets of the vascular domain.
* `-geo_file` specifies the vtu file for the volumetric mesh.
* `-sur_file_in` specifies the vtp or vtu file for the inlet surface.
