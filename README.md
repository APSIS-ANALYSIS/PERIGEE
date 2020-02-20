# PERIGEE
PERIGEE is a nonlinear dynamic finite element / isogeometric analysis code for multiphysics simulations. The code has been developed with the goal of providing a single C++ framework for parallel implementation of different physics problems using different element technologies.

## Table of Contents

- [Install](#Install)
- [References](#References)

## Install
Please follow the following steps to compile PERIGEE.

1. For Sherlock@Stanford users, go to step 2. Before compiling PERIGEE, one has to install several libraries. A quick guide for library installation is [here](docs/install_external_libs.md) and a slightly advanced guide is [there](docs/install-advanced.md). After the libraries are all properly installed, proceed to step 3.

2. On Sherlock@Stanford, all the libraries have been installed in `/home/groups/amarsden/lib-perigee`. You need to add the following to your `.bash_profile` or `.bashrc` file to define `MACHINE_NAME` as an environment variable, and then proceed to step 3. With the `MACHINE_NAME`, PERIGEE can load the proper configuration file for compiling.
```sh
export MACHINE_NAME=sherlock
export LD_LIBRARY_PATH=/home/groups/amarsden/lib-perigee/VTK-7.1.1-shared/lib:$LD_LIBRARY_PATH
```
 
3. After the libraries are installed, one has to modify the configuration file in the [conf](conf) folder, following the steps [here](docs/configure_perigee_guide.md). *If you are on Sherlock@Stanford, you do not need to do anything at this step. The CMake configuration file for Sherlock is [here](conf/stanford_sherlock.cmake). As long as you have your machine named as `sherlock`, PERIGEE will load the proper CMake file and compile the code*.

## Build
We need to first create a folder `build` out of the PERGIEE source directory. Enter the folder, and run the following command to build a suite of heat equation solvers.
```sh
CMake ~/PERIGEE/examples/nonlinear_heat_eqn/
make
```

## References
### Finite Element Method
* [The Finite Element Method: Linear Static and Dynamic Finite Element Analysis](https://www.amazon.com/Finite-Element-Method-Mechanical-Engineering/dp/0486411818/ref=sr_1_2?keywords=the+finite+element+method&qid=1566093145&s=books&sr=1-2) by Thomas J.R. Hughes
* Incompressible Flow and the Finite Element Method, Volume 1: Advection-Diffusion and Isothermal Laminar Flow by P.M. Gresho and R.L. Sani

### C++
* [www.learncpp.com](http://www.learncpp.com).
* [Google C++ Style](https://google.github.io/styleguide/cppguide.html).

## Contact
Ju Liu, liujuy@gmail.com, liuju@stanford.edu
