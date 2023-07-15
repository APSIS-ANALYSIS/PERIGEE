# Instruction
This folder contains a suite of codes for performing patient-specific vascular FSI simulations using a small-strain membrane wall model.

## Table of Contents

- [Preprocessor](#Preprocessor)
- [FSI analysis](#FSI-Analysis)
- [Wall prestress generation](#Wall-prestress-generation)
- [Simulation pipeline](#Simulation-pipeline)
- [Postprocessing](#Postprocessing)

## Preprocessor
In the [preprocessor](preprocess_sv_tets.cpp), one needs to prepare the mesh, assign proper boundary conditions, and do mesh partitioning. The following arguments determine the basic setting of the geometry as well as the mesh partitioning.
* `-cpu_size` determines the number of CPUs you want to run for the analysis solver.
* `-elem_type` determines the type of mesh you want to use.
* `-num_outlet` determines the number of outlets of the vascular domain.
* `-geo_file` specifies the vtu file for the volumetric mesh.
* `-sur_file_in` specifies the vtp or vtu file for the inlet surface.
* `-sur_file_wall` specifies the vtp or vtu file for the wall surface.
* `-sur_file_out_base` specifies the vtp or vtu file for the outlet surfaces. Together with the number of outlets, the code will load sur_file_out_base + xxx.vtp (vtu) from the disk.

The following arguments determine the wall thickness and Young's modulus.
* `-is_uniform_wall` is a bool argument that determines if we want to have uniform wall properties. It determines which constructor is called for ElemBC_3D_tet_wall. If it is false, the wall properties are determined from a centerline file. The file for the centerline should be named as `centerlines.vtp` and be placed in the preprocessor job folder. The wall thickness is 20 percent of the radius. The Young's modulus is calculated from an empirical formula, which is. defined in ElemBC_3D_tet_wall::compute_youngsmod. If the argument's value is true, use `-wall_thickness` and `-wall_youngsmod` to define the wall thickness and Young's modulus, respectively.
* `-wall_springconst` and `-wall_dampingconst` define the supporting tissue model's elastic coefficient and damping coefficient, respectively.

We recommend using `-is_uniform_wall NO`, so that one can generate more physiologically realistic wall properties. Note: make sure that the correct centerlines.vtp is placed in the job folder.

So far, the thickness-to-radius ratio, the formula for the Young's modulus, the supporting tissue's two parameters are assumed to take a uniform value. In a more sophsiticated manner, one may set them with different values in different regions of the wall. To do so, one needs to go over the thrid constructor for the ElemBC_3D_tet_wall class.

The constructor of ElemBC_3D_tet_wall class will always write `varwallprop.vtp` as a wall surface file, with wall properties as cell or pint attributes. 

The following argument determines the boundary condition type on the inlet and outlets.
* `-cmmbc_type` is an integer flag that determines NodalBC_3D_CMM. If it is `0` (default value), the wall is set to be deformable; if it is `1`, the wall is set as a homogeneous boundary meaning the wall is rigid; if it is `2`, the variables in the fluid subdomain is fixed so that one can solve the wall mechanics for prestress.
* `-ringbc_type` is an integer flag that determines NodalBC_3D_ring, which also affects NodalBC_3D_CMM. If it is `0` (default value), the ring nodes are fully clamped; if it is `1`, the ring nodes are allowed to move within their original plane. The ringbc_type only affects the boundary condition when the cmmbc_type = 0 or 2.

## FSI Analysis
In the actual [analysis driver](cmm_driver.cpp), one can perform both deformable membrane FSI analysis (i.e., CMM-FSI) and the rigid wall CFD analysis. There are a few things that we want to mention about its arguments.

* `-fl_density` and `-fl-mu` define the fluid density and viscosity, respectively.
* `-wall_density`, `-wall_poisson`, and `-wall_kappa` define the wall density, Poisson's ratio, and the shear correction factor. We assume they are uniform.
* `-nqp_tet` and `-nqp_tri` defines the number of volume and surface quadrature points. Note that if you use quadratic mesh, be sure to set them to be `29` and `13`, respectively, otherwise the code will be killed.
* `-inflow_file` specifies the inflow in terms of fourier series, with the default value `inflow_fourier_series.txt`. If this file is placed in the job folder and has correct format, the code will load it to generate a (pulsatile) inflow profile. Otherwise, the code will set the initial flow rate to be zero and increase the flow rate with respect to time until reaching a target flow rate value. The linear incremental rate and the target flow rate can be controlled by `-inflow_thd_time` and `-inflow_tgt_rate`.
* `-lpn_file` defines the LPN data for the outlets. Therefore, one should place the file in the job folder, and the default name for the file is `lpn_rcr_input.txt`.
* `-nz_estimate` defines the estimation of num nonzeros per row for the sparse tangent matrix.
* `-bs_beta` defines backflow stabiization parameter.
* `-rho_inf` defines generalized-alpha parameter.
* `-is_backward_Euler` defines whether adopts the backward_Euler method.
* `-c_tauc` defines scaling factor for tau_c：0.0, 0.125, or 1.0.
* `-inflow_thd_time` defines the time for linearly increasing inflow to reach steady state.
* `-inflow_tgt_rate` defines the rate for linearly increasing inflow to reach steady state.
* `-part_file` defines the base mesh partition files from preprocess program.
* `-nl_rtol` defines the convergence criterion relative tolerance.
* `nl_atol` defines the convergence criterion absolute tolerance.
* `-nl_maxits` defines the maximum number of nonlinear iterations
* `-nl_refreq` defines the frequency of tangent matrix renewal.
* `-nl_threshold` defines the threshold of tangent matrix renewal.
* `-init_time` defines the time of initial condition.
* `-fina_time` defines the end time of simulation.
* `-init_step` defines the time step.
* `-init_index` defines the index of initial condition.
* `-ttan_freq` defines the frequency of tangent matrix renewal.
* `-sol_rec_freq` defines the frequency for recording the solution.
* `-sol_name` defines the base name of the solution file.
* `-is_restart` defines whether run restart process. Restart here means the effect of wall deformation is imposed on the fluid simulation.
* `-restart_index` defines the restart solution time index.
* `-restart_time` defines the restart time.
* `-restart_step` deifnes the restart simulation time step size.
* `-restart_name` defines the restart solution base name.

Make sure that the HDF5 files, the inflow file, and the LPN file are placed in the analysis job folder.

## Wall prestress generation
The [wall solver](wall_solver.cpp) generate the prestress in the arterial wall. This solver requires the following files.
* The solver command line input `solver_cmd.h5`, which will provide the material properties set in the CFD analysis.
* The mesh partitioning file `part_xxxxx.h5`, generate by the preprocessor with `-cmmbc_type 1`.
* The CFD steady-state solution saved as `SOL_re`.

Therefore, the user should be copy the `solver_cmd.h5` and `SOL_re` from the CFD analysis and rerun the preprocessor to generate the `part_xxxxx.h5` files for prestress generation. 

A few things about its input arguments:
* `-prestress_disp_tol` determines the stoping criterion for the fixed-point iteration, with the default value being `1.0e-6`.
* `-init_step` sets the time step size. As we are driving for a steady-state solution, it is recommended to make the time step large.
* `-is_backward_Euler` tells this solver if we want to use Backward Euler or the Generalized-alpha method. It is recommended to use the Backward Euler method.
* `-is_record_sol` tells the wall solver if we want to save the solution files.

## Simulation pipeline
Here we describe a standard pipeline for performing CMM-FSI simulations. 

1. One may want to prepare the mesh partitioning for a rigid wall CFD simulation. We assume that the mesh files whole_vol.vtu, inflow_vol.vtp, wall_vol.vtp, outflow_vol_xxx, and the centerline file centerlines.vtp have been placed in the job folder.
```sh
./preprocess3d -cmmbc_type 1 -ringbc_type 0 -is_uniform_wall NO -num_outlet 46 -cpu_size 60 -elem_type 501
```
2. With the generated HDF5 files as well as the LPN file `lpn_rcr_input.txt`, one may call the analysis solver to generate a flow profile at the diastolic phase.
```sh
mpirun -np 60 ./cmm_tet_3d \
  -fl_density 1.00 -fl_mu 4.0e-2 -wall_density 1.0 -wall_poisson 0.5 \
  -nqp_tet 5 -nqp_tri 4 -init_step 4.055e-3 -fina_time 4.866 -is_backward_Euler YES \
  -inflow_file inflow_fourier_series_steady.txt -lpn_file lpn_rcr_input.txt \
  -nl_refreq 1 -nl_rtol 1.0e-6 -nl_atol 1.0e-6 -nl_dtol 1.0e8 -nl_maxits 20 \
  -ttan_freq 100 -sol_rec_freq 10 \
  -is_restart NO \
  -ksp_type gmres -pc_type asm -ksp_rtol 1.0e-2 -ksp_atol 1.0e-50 -ksp_max_it 200 -ksp_gmres_restart 200 \
  -log_view
```
   * In the above script, one prepares a file `inflow_fourier_series_steady.txt` that has only one frequency with coefficient `a_0` matches the diastolic inflow.
   
   * Otherwise, without preparing the above text file, one could also set `-inflow_thd_time` be smaller than the `-fina_time` and `-inflow_tgt_rate` equaling the diastolic inflow rate. The flow will be linearly ramped up to that value.

3. If the flow solver runs successfully, there will be a steady-state solution recorded and we rename the solution files (for example `SOL_900001000`, `dot_SOL_900001000`, `SOL_disp_900001000`, and `dot_SOL_disp_900001000` if the last step index is 1000) to `SOL_re`, `dot_SOL_re`, `SOL_disp_re`, and `dot_SOL_disp_re`. Copy `SOL_re` and `solver_cmd.h5` to the folder for the prestress generation. Run the following to generate the `prestress_pxxxxx.h5` files.
```sh
./preprocess3d -cmmbc_type 2 -ringbc_type 0 -is_uniform_wall NO -num_outlet 46 -cpu_size 60 -elem_type 501

mpirun -np 60 ./wall_solver \
  -is_record_sol NO -is_backward_Euler YES -prestress_disp_tol 1.0e-6 \
  -init_step 1.0e-2 -fina_time 1.0 \
  -nl_refreq 1 -nl_rtol 1.0e-6 -nl_atol 1.0e-6 -nl_dtol 1.0e8 -nl_maxits 20 \
  -ttan_freq 100 -sol_rec_freq 1 \
  -ksp_type gmres -pc_type asm -ksp_rtol 1.0e-2 -ksp_atol 1.0e-50 -ksp_max_it 200 -ksp_gmres_restart 200 \
  -log_view
```
Note, in the generation of prestress, it is recommended to clamp the ring (-ringbc_type 0), and not prescribing any supporting tissue boundary condition. This can be critical for quick convergence of the prestress solver.

4. Now we can rerun the preprocessor to assign boundary conditions for the deformable wall simulation.
```sh
./preprocess3d -cmmbc_type 0 -ringbc_type 0 -is_uniform_wall NO -num_outlet 46 -cpu_size 60 -elem_type 501
```

5. Now make sure that the HDF5 files generated from the above script, the prestress_pxxxxx.h5 file, the LPN file `lpn_rcr_input.txt`, and the inflow file `inflow_fourier_series.txt` have all been placed in the same job folder. Also, one needs to copy all the steady-state solutions from the previous CFD analysis and place them with the SOL_re file in the job folder. Run the following for the CMM simulation.
```sh
mpirun -np 60 ./cmm_tet_3d \
   -nz_estimate 1000 -fl_density 1.00 -fl_mu 4.0e-2 -wall_density 1.0 -wall_poisson 0.5 \
   -nqp_tet 5 -nqp_tri 4 \
   -fina_time 2.433 \
   -lpn_file lpn_rcr_input.txt \
   -nl_refreq 1 -nl_rtol 1.0e-6 -nl_atol 1.0e-6 -nl_dtol 1.0e8 -nl_maxits 20 \
   -ttan_freq 100 -sol_rec_freq 40 \
   -is_restart YES \
   -restart_index 0 \
   -restart_time 0.0 \
   -restart_step 4.055e-4 \
   -restart_name SOL_re \
   -restart_disp_name SOL_disp_re \
   -ksp_type gmres -pc_type asm -ksp_rtol 1.0e-2 -ksp_atol 1.0e-50 -ksp_max_it 200 -ksp_gmres_restart 200 \
   -log_view
```
Notice that we want to run this as a restart run with the restart solution coming from the steady state CFD analysis, since the wall loads the prestress, which is  balancing with the fluid pressure.

## Postprocessing
To run postprocessing scripts, one has to repartition the mesh. Running this script is fairly similar to the preprocessing. One just needs to specify the number of CPUs for postprocessing. Make sure that the original geometry files (in vtu or vtp format) as well as the `preprocessor_cmd.h5` are in the job folder, and then run the following, which means repartition the mesh into 16 subdomains.
```sh
./prepost -cpu_size 16
```
Next, as an example of postprocessing, we visualize the data. First, we need to make sure that we have all the necessary data, including `epart.h5`, `node_mapping.h5`, all postprocessing mesh partitioning generated hdf5 files, and solution files (e.g. SOL_900000000). Then run the following to do the visualization.
```sh
mpirun -np 16 ./vis_cmm -time_start 1000 -time_step 1000 -time_end 5000
```
The above command convert the solution files to vtk files, from time step 1000 to 5000, with incremental step being 1000.
