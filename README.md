# PERIGEE
PERIGEE is a nonlinear dynamic finite element analysis code for multiphysics simulations. The code has been developed with the goal of providing an object-oriented framework for parallel implementation of multiphysics problems. Copyright and licensing information can be found in files [LICENSE](LICENSE).

## Table of Contents

- [Install](#Install)
- [Build](#Build)
- [Tutorial](#Tutorial)
- [Simulation Samples](#Simulation-Samples)
- [References](#References)

## Install
Please go through the following steps to install external libraries.

1. Before compiling PERIGEE, you will have to install a few external libraries. A quick guide for library installation is [here](docs/install_external_libs.md) and a more advanced guide is [there](docs/install-advanced.md). After the libraries are all properly installed, proceed to step 3.

2. You need to add the following to your `.bash_profile` or `.bashrc` file to define `MACHINE_NAME` as an environment variable, and then proceed to step 3. For example, I named one of my machine as `sherlock`, and PERIGEE will detech the `MACHINE_NAME` and load an appropriate configuration file.
```sh
export MACHINE_NAME=sherlock
```
Additionally, because VTK is typically installed as a shared library in a non-standard folder, one has to edit the `LD_LIBRARY_PATH` environmental variable for the linker. So you will have to edit the `LD_LIBRARY_PATH` variable as above to tell the computer how to locate the VTK libraries. For more information on this environmental variable, see [here](http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html).

```sh
export LD_LIBRARY_PATH=/home/groups/amarsden/lib-perigee/VTK-7.1.1-shared/lib:$LD_LIBRARY_PATH
```
 
3. After the libraries are installed, add or modify the configuration file in the [conf](conf) folder, following the steps [here](docs/configure_perigee_guide.md).


## Build
First, create a folder `build` out of the source directory. Enter that folder, and run the following commands to build, as an example, a suite of heat equation solvers.
```sh
CMake ~/PERIGEE/examples/nonlinear_heat_eqn/
```
CMake will print some information on screen. Pay a look at the variable `CMAKE_BUILD_TYPE`. If its value is `Debug`, this means your code will be compiled in the debug mode. If you are not developing the code and wants to make the code faster, run CMake as follows,
```sh
CMake ~/PERIGEE/examples/nonlinear_heat_eqn/ -DCMAKE_BUILD_TYPE=Release
```
Now the value of `CMAKE_BUILD_TYPE` is set to `Release`. The code will be compiled in the optimized mode. For more information about the compiler, please refer to [this](https://stackoverflow.com/questions/48754619/what-are-cmake-build-type-debug-release-relwithdebinfo-and-minsizerel/48755129). Of course, a fully optimized code requires that your external libraries, especially PETSc, are compiled in the optimized mode. See the [advanced guide](docs/install-advanced.md) to learn how to build libraries in optimal mode. CMake will generate the Makefile for you and you just need to run the following command to compile the source code.
```sh
make
```
Of course you may add `-j2` to run Makefile with 2 threads. If the make complains about the auto keyword or the nullptr, your default compiler does not support C++11. You may add `SET(CMAKE_CXX_STANDARD 11)` in your .cmake configuration file to enforce the C++11 standard. 

## Tutorial
In general, one has to go through the following steps for simulation.
* Obtain the mesh in vtu/vtp format from SimVascular or Gmsh.
* Run a preprocessor to load the mesh, assign boundary conditions, and partition the mesh. The preprocessor is a *serial* code and may need to be run on a large memory cluster node if you are dealing with a very large problem.
* Run a finite element analysis code to solve the partial differential equations. The solutions will be saved on disk in the binary format.
* Run a preprocessor for postprocessing. This step re-partition the mesh to make preparations for postprocessing, such as visualization, error calculation, etc. Similar to the preprocessor, this routine should be run in *serial* and may consume a lot memory if your mesh is fine. With this routine, we are able to run the postprocessing routine with different number of CPUs. For example, we run FEM analysis with, say, 360 CPUs; visualizing the solution is much less intensive in computing and may only need, say, 24 CPUs. So you should repartition the domain into 24 sub-domains in this step.
* Run a postprocessor in parallel. Often, this step refers to the visualization of the solutions. The visualzation routine will read the binary solution files and write the data into (parallel) vtu/vtp format. Then the data can be visualized in Paraview.

## Simulation Samples
The pulmonary model is built from a healthy 20 month old male, consisting of 772 outlets. A rigid wall CFD simulation is performed with 26 million linear elements. The model is prepared by M. Dong and the mesh is generated by W. Yang. Click the image for a youtube video.

Reference: J. Liu, W. Yang, M. Dong, and A.L. Marsden, "The nested block preconditioning technique for the incompressible Navier-Stokes equations with emphasis on hemodynamic simulations.", Computer Methods in Applied Mechanics and Engineering, 2020.

[![Pulmonary CFD](http://img.youtube.com/vi/nbrhpyRE4IU/0.jpg)](https://www.youtube.com/watch?v=nbrhpyRE4IU "Pulmonary CFD")

A fluid-structure interaction simulation of a pulmonary model is performed using **the unified continuum and variational multiscale formulation**. The model and mesh are prepared by W. Yang. The solid model is *fully incompressible* and is numerically modeled via the residual-based variational multiscale formulation.

Reference: J. Liu, W. Yang, I.S. Lan, and A.L. Marsden, "Fluid-structure interaction modeling of blood flow in the pulmonary arteries using the unified continuum and variational multiscale formulation.", Mechanics Research Communications, 2020.

[![Pulmonary FSI](http://img.youtube.com/vi/Y84vSN64ZCk/0.jpg)](https://www.youtube.com/watch?v=Y84vSN64ZCk "Pulmonary FSI")

## References
### Theory
* J. Liu and A.L. Marsden, [A unified continuum and variational multiscale formulation for fluids, solids, and fluid-structure interaction](https://doi.org/10.1016/j.cma.2018.03.045), Computer Methods in Applied Mechanics and Engineering, 337:549-597, 2018.
* I.S. Lan, J. Liu, W. Yang, and A.L. Marsden, [A reduced unified continuum formulation for vascular fluid-structure interaction](https://doi.org/10.1016/j.cma.2022.114852), Computer Methods in Applied Mechanics and Engineering, 394:114852, 2022.

### HPC
* D. Goldberg, What every computer scientist should know about floating-point arithmetic.
* U. Drepper, What every programmer should know about memory.

### C++
* [www.learncpp.com](http://www.learncpp.com).
* [Google C++ Style](https://google.github.io/styleguide/cppguide.html).

## Contact
Ju Liu, liujuy@gmail.com, liuj36@sustech.edu.cn
