# Instruction
This folder contains a suite of codes for performing patient-specific vascular FSI simulations using a small-strain membrane wall model.

## Table of Contents

- [Preprocessor](#Preprocessor)
- [FSI analysis](#FSI-Analysis)

## Preprocessor
In the [preprocessor](preprocess_sv_tets.cpp), one needs to prepare the mesh, assign proper boundary conditions, and do mesh partitioning. There are a few critical input arugments. The following arguments determine the basic setting of the geometry as well as the mesh partitioning.
* `-cpu_size` determines the number of CPUs you want to run for the analysis solver.
* `-elem_type` determines the type of mesh you want to use.
* `-num_outlet` determines the number of outlets of the vascular domain.
* `-geo_file` specifies the vtu file for the volumetric mesh.
* `-sur_file_in` specifies the vtp or vtu file for the inlet surface.
* `-sur_file_wall` specifies the vtp or vtu file for the wall surface.
* `-sur_file_out_base` specifies the vtp or vtu file for the outlet surfaces. Together with the number of outlets, the code will load sur_file_out_base + xxx.vtp (vtu) from the disk.

The following arguments determine the wall thickness and Young's modulus.
* `-is_uniform_wall` is a bool argument that determines if we want to have uniform wall properties. It determines which constructor is called for ElemBC_3D_tet_wall. If it is false, the wall properties are determined from a centerline file. The file for the centerline should be named as centerlines.vtp and be placed in the preprocessor folder. The wall thickness is 20 percent of the radius. The Young's modulus is calculated from an empirical formula, which is. defined in ElemBC_3D_tet_wall::compute_youngsmod. If it is true, the following two arguments determien the actual properties.
    * `-wall_thickness` gives the wall thickness.
    * `-wall_youngsmod` gives the Young's modulus of the wall.
* `-wall_springconst` defines the supporting tissue model's elastic coefficient.
* `-wall_dampingconst` defines the supporting tissue model's damping coefficient.

We recommend using setting this argument to false, so that one can generate more physiologically realistic wall properties. Also, make sure that the correct centerlines.vtp is placed in the same job folder.

So far, the thickness-to-radius ratio, the formula for the Young's modulus, the supporting tissue's two parameters are assumed to take a uniform value. In a more sophsiticated manner, one may set them with different values in different regions of the wall. To do so, one needs to go over the thrid constructor for the ElemBC_3D_tet_wall class.

The following argument determines the boundary conditions on the inlet and outlets.
* `-cmmbc_type` is an integer flag that determines NodalBC_3D_CMM. If it is `0` (default value), the wall is set to be deformable; if it is `1`, the wall is set as a homogeneous boundary meaning the wall is rigid; if it is `2`, the variables in the fluid subdomain is fixed so that one can solve the wall mechanics for prestress.
* `-ringbc_type` is an integer flag that determines NodalBC_3D_ring, which also affects NodalBC_3D_CMM. If it is `0` (default value), the ring nodes are fully clamped; if it is `1`, the ring nodes are allowed to move within their original plane.

## FSI Analysis
In the actual [analysis driver](cmm_driver.cpp), one can perform both deformable membrane type FSI analysis (i.e., CMM-FSI) and the rigid wall CFD analysis. There are a few things that we want to mention about its arguments.

First, the actual fluid and solid material properties used for the simulation are defined through the following arguments.
* `-fl_density` and `-fl-mu` defines the fluid density and viscosity.
* `-wall_density`, `-wall_poisson`, and `-wall_kappa` defines the wall density, Poisson's ratio, and the shear correction factor. We assume they are uniform in the wall.
