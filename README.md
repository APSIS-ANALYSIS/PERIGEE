<img src="./docs/PERIGEE-logo.png"  width="36%" height="36%"> 

PERIGEE is a nonlinear dynamic finite element analysis code for multiphysics analysis. The code has been developed with the goal of providing an object-oriented framework for parallel implementation of multiphysics problems. Copyright and licensing information can be found in [LICENSE](LICENSE).

## Table of Contents

- [Install](#Install)
- [Build](#Build)
- [Tutorial](#Tutorial)
- [Contributing](#Contributing)
- [Simulation Samples](#Simulation-Samples)
- [References](#References)

## Install
We recommend using UNIX-like operating systems, such as Linux or MacOS, for the code development. If you are a Windows user, you may refer to [this](docs/install_guidance_windows.md) for a detailed install guide. The following instructions are based on a Linux Ubuntu system, and there could be minor differences for Mac systems.

1. A quick guide for library installation is [here](docs/install_external_libs.md) and a more advanced guide is [there](docs/install-advanced.md). After the libraries are all properly installed, proceed to step 2.

Notice that VTK is typically installed as a shared library in a non-standard folder. One therefore has to edit the `LD_LIBRARY_PATH` environmental variable for the linker to locate the .so files. Open the `.bash_profile` or `.bashrc` file and edit the `LD_LIBRARY_PATH` variable. See below is an example with my VTK installed at `/Users/juliu/lib/VTK-8.2.0/`.

```sh
export LD_LIBRARY_PATH=/Users/juliu/lib/VTK-8.2.0/lib:$LD_LIBRARY_PATH
```
For more information on this environmental variable, see [here](http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html).

2. After the libraries are installed, add a configuration file named as `system_lib_loading.cmake` in the [conf](conf) folder. You may find a file called `system_lib_loading_example.cmake`, which is an example. In this file, you will have to specify the paths for the external libraries,

 - Set `VTK_DIR` to the VTK library location (e.g. `/home/jliu/lib/VTK-7.1.1-shared`).
 - Set `PETSC_DIR` to the PETSc library location (e.g. `/home/jliu/lib/petsc-3.11.3`).
 - Set `PETSC_ARCH` to the value used in PETSc installation (e.g. `arch-linux2-c-debug`).
 - Set `METIS_DIR` to the METIS library location (e.g. `/home/jliu/lib/metis-5.0.3`).
 - Set `HDF5_DIR` to the HDF5 library location (e.g. `/home/jliu/lib/hdf5-1.8.16`).
 - Set `CMAKE_C_COMPILER` to `$PETSC_DIR/$PETSC_ARCH/bin/mpicc`
 - Set `CMAKE_CXX_COMPILER` to `$PETSC_DIR/$PETSC_ARCH/bin/mpicxx`

After the edit, save the CMake file and rename it as `system_lib_loading.cmake`, and you have your own configuration file set up. Notice that we have the file name `system_lib_loading.cmake` added in .gitignore, meaning that git will not track this file. You may want to keep a copy of this file out of PERIGEE, because when you switch to other branches, PERIGEE may not keep a copy of this file. 

## Build
First, create a folder `build` out of the source directory. Enter that folder, and run the following commands to build, as an example, a suite of heat equation solvers.
```sh
CMake ~/PERIGEE/examples/linearPDE/
```
CMake will print some information on the screen. Pay a look at the variable `CMAKE_BUILD_TYPE`. If its value is `Debug`, this means your code will be compiled in the debug mode. If you want to make the code faster, run CMake as follows,
```sh
CMake ~/PERIGEE/examples/linearPDE/ -DCMAKE_BUILD_TYPE=Release
```
Now the value of `CMAKE_BUILD_TYPE` is set to `Release`. The code will be compiled in the optimized mode. For more information about the compiler, please refer to [this](https://stackoverflow.com/questions/48754619/what-are-cmake-build-type-debug-release-relwithdebinfo-and-minsizerel/48755129). Of course, a fully optimized code requires that your external libraries, especially PETSc, are compiled in the optimized mode also. Refer to the [advanced guide](docs/install-advanced.md) for more info on building libraries in a release mode. After CMake generates the Makefile for you, you need to run the following command to compile the source code.
```sh
make
```
Of course, you may add `-j2` to run Makefile with 2 threads. If the make complains about the auto keyword or the nullptr, your default compiler does not support C++11. You may add `SET(CMAKE_CXX_STANDARD 11)` in your .cmake configuration file to enforce the C++11 standard. 

## Tutorial
In general, one has to go through the following steps for simulation.
* Obtain the mesh in vtu/vtp format from a front-end code, e.g., SimVascular or Gmsh.
* Run a preprocessor to load the mesh, assign boundary conditions, and partition the mesh. The preprocessor is a *serial* code and may need to be run on a large memory cluster node if you are dealing with a very large problem.
* Run a finite element analysis code to solve the partial differential equations. The solutions will be saved on disk in the binary format.
* Run a preprocessor for postprocessing. This step re-partition the mesh to make preparations for postprocessing, such as visualization, error calculation, etc. Similar to the preprocessor, this routine should be run in *serial* and may consume a lot memory if your mesh is fine. With this routine, we are able to run the postprocessing routine with different number of CPUs. For example, we run FEM analysis with, say, 360 CPUs; visualizing the solution is much less intensive in computing and may only need, say, 24 CPUs. So you should repartition the domain into 24 sub-domains in this step.
* Run a postprocessor in parallel. Often, this step refers to the visualization of the solutions. The visualzation routine will read the binary solution files and write the data into (parallel) vtu/vtp format. Then the data can be visualized in Paraview.

## Contributing
We welcome contributions from the community! If you are interested in contributing to this project, please follow the guidelines below.

### Code Style
Please follow these guidelines when contributing to the codebase:

- Follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
- Make sure your code is well-documented with comments explaining the purpose and functionality.
- Write clear and concise commit messages.
- Keep your changes focused on a single feature or bug fix.

### Testing
- Please make sure that your changes are covered by tests.
- Run all tests to verify that your changes do not break existing functionality.
- If adding new features, add appropriate tests to cover the new functionality.

### Issues and Bugs
- If you find a bug, please report it by opening an issue.
- Provide as much detail as possible, including steps to reproduce the issue and any relevant error messages or logs.

### Pull Request Process
1. Ensure your changes pass all tests and adhere to the project's coding standards.
2. Update the README.md file if necessary to reflect your changes.
3. Open a pull request and provide a detailed description of your changes.
4. A project maintainer will review your pull request and provide feedback.
5. Once approved, your changes will be merged into the main branch.

### Community
- Be respectful and considerate in your interactions with others.
- Follow the [Code of Conduct](CODE_OF_CONDUCT.md) to ensure a welcoming and inclusive environment for all contributors.

## Simulation Samples
The vortex-induced vibration of an elastic plate with Re $\approx$ 3 $\times$ 10^4. The mesh consists of 18 million linear tetrahedral elements for the fluid and 0.7 million elements for the solid. The variational multiscale formulation provides the LES technique in the flow problem, and the time integration is based on the generalized-&alpha; scheme.

[![Pulmonary CFD](http://img.youtube.com/vi/QiSkyBMGhmI/0.jpg)](https://www.youtube.com/watch?v=QiSkyBMGhmI "Vortex induced vibration")

A fluid-structure interaction simulation of a pulmonary model is performed using **the unified continuum and variational multiscale formulation**. The model and mesh are prepared by W. Yang. The solid model is *fully incompressible* and is numerically modeled via the residual-based variational multiscale formulation.

[![Pulmonary FSI](http://img.youtube.com/vi/Y84vSN64ZCk/0.jpg)](https://www.youtube.com/watch?v=Y84vSN64ZCk "Pulmonary FSI")

## References
### Theory
* J. Liu and A.L. Marsden, [A unified continuum and variational multiscale formulation for fluids, solids, and fluid-structure interaction](https://doi.org/10.1016/j.cma.2018.03.045), Computer Methods in Applied Mechanics and Engineering, 337:549-597, 2018.
* I.S. Lan, J. Liu, W. Yang, and A.L. Marsden, [A reduced unified continuum formulation for vascular fluid-structure interaction](https://doi.org/10.1016/j.cma.2022.114852), Computer Methods in Applied Mechanics and Engineering, 394:114852, 2022.
* J. Liu, I.S. Lan, O.Z. Tikenogullari, and A.L. Marsden, [A note on the accuracy of the generalized-α scheme for the incompressible Navier-Stokes equations](https://doi.org/10.1002/nme.6550), International Journal for Numerical Methods in Engineering, 122:638-651, 2021.

### Solver technology
* J. Liu, W. Yang, M. Dong, and A.L. Marsden, [The nested block preconditioning technique for the incompressible Navier-Stokes equations with emphasis on hemodynamic simulations](https://doi.org/10.1016/j.cma.2020.113122), Computer Methods in Applied Mechanics and Engineering, 367:113122, 2020.
* Y. Sun, Q. Lu, and J. Liu, [A parallel solver framework for fully implicit monolithic fluid-structure interaction](https://doi.org/10.1007/s10409-024-24074-x), Acta Mechanica Sinica 40:324074, 2024. 

### Verification & Validation
* I.S. Lan, J. Liu, W. Yang, J. Zimmermann, D.B. Ennis, and A.L. Marsden, [Validation of the reduced unified continuum formulation against in vitro 4D-flow MRI](https://doi.org/10.1007/s10439-022-03038-4), Annals of Biomedical Engineering, 51:377-393, 2023.
* J. Liu, W. Yang, I.S. Lan, and A.L. Marsden, [Fluid-structure interaction modeling of blood flow in the pulmonary arteries using the unified continuum and variational multiscale formulation](https://doi.org/10.1016/j.mechrescom.2020.103556), Mechanics Research Communications, 107:103556, 2020.

### Application
* Y. Sun, J. Huang, Q. Lu, X. Yue, X. Huang, W. He, Y. Shi, and J. Liu, [Modeling fibrous tissue in vascular fluid-structure interaction: a morphology-based pipeline and biomechanical significance](https://doi.org/10.1002/cnm.3892), International Journal for Numerical Methods in Biomedical Engineering 41(1):e3892, 2025.
* X. Yue, J. Huang, and J. Liu, [Fluid-structure interaction analysis for abdominal aortic aneurysms: the role of multi-layered tissue architecture and intraluminal thrombus](https://doi.org/10.3389/fbioe.2025.1519608), Frontiers in Bioengineering and Biotechnology 13:1519608, 2025. 

### HPC
* D. Goldberg, What every computer scientist should know about floating-point arithmetic.
* U. Drepper, What every programmer should know about memory.

### C++
* [www.learncpp.com](http://www.learncpp.com).
* [Google C++ Style](https://google.github.io/styleguide/cppguide.html).
* [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)

## Contact
[Dr. Ju Liu](https://ju-liu.github.io), liujuy@gmail.com, liuj36@sustech.edu.cn

## Acknowledgement
National Natural Science Foundation of China, Grant numbers 12172160, 12472201

Shenzhen Science and Technology Program, Grant number JCYJ20220818100600002

<img src="./docs/NSFC_logo.png"  width="22%" height="22%"> 
