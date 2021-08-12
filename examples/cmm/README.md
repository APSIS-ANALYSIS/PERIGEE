# Instruction
This folder contains a suite of codes for performing patient-specific vascular FSI simulations using a small-strain membrane wall model.

## Table of Contents

- [Preprocessor](#Preprocessor)
- [FSI analysis](#FSI-Analysis)
- [Wall prestress generation](#Wall-prestress-generation)

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
    * `-wall_thickness` and `-wall_youngsmod` define the wall thickness and Young's modulus, respectively.
* `-wall_springconst` and `-wall_dampingconst` define the supporting tissue model's elastic coefficient and damping coefficient, respectively.

We recommend using setting this argument to false, so that one can generate more physiologically realistic wall properties. Also, make sure that the correct centerlines.vtp is placed in the same job folder.

So far, the thickness-to-radius ratio, the formula for the Young's modulus, the supporting tissue's two parameters are assumed to take a uniform value. In a more sophsiticated manner, one may set them with different values in different regions of the wall. To do so, one needs to go over the thrid constructor for the ElemBC_3D_tet_wall class.

The following argument determines the boundary conditions on the inlet and outlets.
* `-cmmbc_type` is an integer flag that determines NodalBC_3D_CMM. If it is `0` (default value), the wall is set to be deformable; if it is `1`, the wall is set as a homogeneous boundary meaning the wall is rigid; if it is `2`, the variables in the fluid subdomain is fixed so that one can solve the wall mechanics for prestress.
* `-ringbc_type` is an integer flag that determines NodalBC_3D_ring, which also affects NodalBC_3D_CMM. If it is `0` (default value), the ring nodes are fully clamped; if it is `1`, the ring nodes are allowed to move within their original plane. The ringbc_type only affects the boundary condition when the cmmbc_type = 0 or 2.

## FSI Analysis
In the actual [analysis driver](cmm_driver.cpp), one can perform both deformable membrane type FSI analysis (i.e., CMM-FSI) and the rigid wall CFD analysis. There are a few things that we want to mention about its arguments.

* `-fl_density` and `-fl-mu` defines the fluid density and viscosity. Notice that the fluid density set in the preprocessor does not affect the analysis solver's definition.
* `-wall_density`, `-wall_poisson`, and `-wall_kappa` defines the wall density, Poisson's ratio, and the shear correction factor. We assume they are uniform in the wall.
* `-nqp_tet` and `-nqp_tri` defines the number of volume and surface quadrature points. Note that if you use quadratic mesh, be sure to set them to be `29` and `13`, respectively, otherwise the code will be killed.
* There are several parameters that define the inflow information. The arugment `-inflow_file` specifies the inflow fourier series, with default value `inflow_fourier_series.txt`. If this file is placed in the job folder and has correct format, the code will load it to generate a (pulsatile) inflow profile. Otherwise, the code will set the initial flow rate to be zero and increase the flow rate with respect to time untill reaching a target flow rate value. The linear incremental rate and the target flow rate can be controlled by `-inflow_thd_time` and `-inflow_tgt_rate`.
* `-lpn_file` defines the LPN data for the outlets. Therefore, one should place the file in the job folder, and the default name for the file is `lpn_rcr_input.txt`.

## Wall prestress generation
In the [wall solver](wall_solver.cpp), one can generate the prestress in the arterial wall. This solver requires the following files.
* The solver command line input `solver_cmd.h5`, which will provide the material properties set in the CFD analysis.
* The mesh partitioning file `part_xxxxx.h5`, generate by the preprocessor with `-cmmbc_type 1`.
* The CFD steady state solution saved as `SOL_re`.

Therefore, the user should be copy the `solver_cmd.h5` and `SOL_re` from the CFD analysis and rerun the preprocessor to generate the `part_xxxxx.h5` files for prestress generation. 

We want to mention a few things about its input arguments.
* The argument `-prestress_disp_tol` determines the stoping criterion for the fixed-point iteration, with default value being `1.0e-6`.
* The argument `-init_step` sets the time step size. As we are driving for a steady state solution, it is recommended to make the time step large.
* The argument `-is_backward_Euler` tells this solver if we want to use Backward Euler or the Generalized-alpha method. It is recommended to use Backward Euler.

