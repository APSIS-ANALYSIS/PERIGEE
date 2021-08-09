# Instruction
This folder contains a suite of codes for performing patient-specific vascular FSI simulations using a small-strain membrane wall model.

## Pre-processing
In the preprocessor (source code is preprocess_sv_tets.cpp), one needs to prepare the mesh, assign proper boundary conditions, and do mesh partitioning. There are a few critical input arugments. The following basic arguments decide the basic setting of the problem mesh and partitioning.
* `-cpu_size` determines the number of CPUs you want to run for the analysis solver.
* `-elem_type` determines the type of mesh you want to use.
* `-num_outlet` determines the number of outlets of the vascular domain.
* `-geo_file` specifies the vtu file for the volumetric mesh.
* `-sur_file_in` specifies the vtp or vtu file for the inlet surface.
* `-sur_file_wall` specifies the vtp or vtu file for the wall surface.
* `-sur_file_out_base` specifies the vtp or vtu file for the outlet surfaces. Together with the number of outlets, the code will load sur_file_out_base + xxx.vtp (vtu) from the disk.
The following arguments determines the wall properties.
* `-is_uniform_wall` is a bool argument that determines if we want to have uniform wall properties. If it is true, the following two arguments determien the actual properties.
* - `-wall_thickness` gives the wall thickness.
* - `-wall_youngsmod` gives the Young's modulus of the wall.
