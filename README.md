# PERIGEE
PERIGEE is a nonlinear dynamic finite element / isogeometric analysis code for multiphysics simulations. The code has been developed with the goal of providing an object-oriented framework for parallel implementation of multiphysics problems using different finite element technologies.

## Table of Contents

- [Install](#Install)
- [Build](#Build)
- [References](#References)

## Install
Please follow the following steps to compile PERIGEE.

1. For Sherlock@Stanford users, directly go to step 2. Before compiling PERIGEE, one has to install several libraries. A quick guide for library installation is [here](docs/install_external_libs.md) and a slightly advanced guide is [there](docs/install-advanced.md). After the libraries are all properly installed, proceed to step 3.

2. On Sherlock@Stanford, all the libraries have been installed in `/home/groups/amarsden/lib-perigee`. You need to add the following to your `.bash_profile` or `.bashrc` file to define `MACHINE_NAME` as an environment variable, and then proceed to step 3. With the `MACHINE_NAME`, PERIGEE can load the proper configuration file for compiling. Because VTK is typically installed as a shared library in a non-standard folder, one has to edit the `LD_LIBRARY_PATH` environmental variable for the linker. For more information on this environmental variable, see [here](http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html).
```sh
export MACHINE_NAME=sherlock
export LD_LIBRARY_PATH=/home/groups/amarsden/lib-perigee/VTK-7.1.1-shared/lib:$LD_LIBRARY_PATH
```
 
3. After the libraries are installed, modify the configuration file in the [conf](conf) folder, following the steps [here](docs/configure_perigee_guide.md). *If you are on Sherlock@Stanford, you do not need to do anything at this step. The CMake configuration file for Sherlock is [here](conf/stanford_sherlock.cmake). As long as you have your machine named as `sherlock`, PERIGEE will load the proper CMake file and compile the code*.

## Build
First, create a folder `build` out of the PERGIEE source directory. Enter that folder, and run the following commands to build, as an example, a suite of heat equation solvers.
```sh
CMake ~/PERIGEE/examples/nonlinear_heat_eqn/
```
CMake will print some information on screen. Pay a look at the variable `CMAKE_BUILD_TYPE`. If its value is `Debug`, this means your code will be compiled in the debug mode. If you want to have better code efficiency, run CMake as follows,
```sh
CMake ~/PERIGEE/examples/nonlinear_heat_eqn/ -DCMAKE_BUILD_TYPE=Release
```
Now the value of `CMAKE_BUILD_TYPE` is set to `Release`. The code will be compiled in the optimized mode. For more information about the compiler, please refer to this [note](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html). CMake will generate the Makefile for you and you just need to run the following command to compile the source code.
```sh
make
```
Of course you may add `-j2` to run Makefile with 2 threads. If the make complains about the auto keyword or the nullptr, your default compiler does not support C++11. You may add `set (CMAKE_CXX_STANDARD 11)` in your .cmake configuration file to enforce the C++11 standard. 

## References
### Finite Element Method
* [The Finite Element Method: Linear Static and Dynamic Finite Element Analysis](https://www.amazon.com/Finite-Element-Method-Mechanical-Engineering/dp/0486411818/ref=sr_1_2?keywords=the+finite+element+method&qid=1566093145&s=books&sr=1-2) by Thomas J.R. Hughes
* Incompressible Flow and the Finite Element Method, Volume 1: Advection-Diffusion and Isothermal Laminar Flow by P.M. Gresho and R.L. Sani

### C++
* [www.learncpp.com](http://www.learncpp.com).
* [Google C++ Style](https://google.github.io/styleguide/cppguide.html).

## Contact
Ju Liu, liujuy@gmail.com, liuju@stanford.edu
