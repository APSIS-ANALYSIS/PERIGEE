# Advanced Installation Guide

The installation of external libraries may involve more settings during the configuration stage. In particular, the installation setting of VTK and PETSc may have a major impact on the overall code performance. Thereby, I provide a detailed external library installation for advanced users.

## Advanced Guide for VTK Installation
The official installation guide for VTK can be found [here](https://vtk.org/Wiki/VTK/Configure_and_Build). The configuration of VTK is done by CMake. Instead of running one-pass cmake command as we given in the [introductary guide](install_external_libs.md), one may first run plain cmake to generate an initial configuration file named CMakeCache.txt.
```sh
$ cmake $HOME/VTK-7.1.1-src
```
Open the `CMakeCache.txt file`, which contains the configuration details of the default build. Pay attention to the following variables.
* `BUILD_SHARED_LIBS` Its default value is `ON`, and instructs the VTK library to be built as shared libraries. The shared libraries ends with .so. If the library is of the shared type, the Linux loder needs to search the library during the run time, and you are responsible of telling the Linux system where the shared libraries are. This is done through defining the environment variable `LD_LIBRARY_PATH`. It is a colon-separated list of paths just like `PATH` in your .bashrc file. If you trun the value of `BUILD_SHARED_LIBS` to `OFF`, the builder will make a static VTK library, and you do not need to worry about `LD_LIBRARY_PATH` in this case.

* `CMAKE_BUILD_TYPE` Its default value is `Debug`. Debug mode is slow but may bring detailed error messages when there are bugs. It is recommended on your local machine where you do the code development. If you turn its value to `Release`, the library will be compiled with optimization flags on. It is recommended for building on clusters.

* `CMAKE_INSTALL_PREFIX` Its default value is `/usr/local`. I do not recommend installing your library in a default system folder. In may cases, you do not have the permission to write in the `/usr/local` folder. Change its value to a place where you have full permissions. In our [guide](install_external_libs.md), we set its value to be `$HOME/lib/VTK-7.1.1-shared`.

After you finished modifying the CMakeCache.txt file, save and close it. Then rerun the cmake command one more times and proceed to run make.
```sh
$ cmake $HOME/VTK-7.1.1-src
$ make -j 6
$ make install
```
## Advanced Guide for PETSc Installation
The PETSc package has an official installation guide [page](https://www.mcs.anl.gov/petsc/documentation/installation.html). The package installation is controlled through the `configure` command.
* `--with-mpi-dir=/home/jliu/lib/mpich-3.2` The `--with-mpi-dir` tells the PETSc that there exists a MPI library installed in the computer, and the PETSc will not have to download and compile a MPICH during installation. My own mpich is installed in `/home/jliu/lib/mpich-3.2`.

* `--with-hypre=1 --download-hypre` This flag tells the PETSc installer to install the [Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) package as well. 

* `--with-debugging=yes` The `--with-debugging` flag is similar to the `CMAKE_BUILD_TYPE` in CMake. It tells the installer if you want to have a debug mode or an optimized mode for the library.

* `--prefix=/home/jliu/lib/petsc-3.11.3-debug` The `--prefix` flag tells the installer where you want to install the PETSc package. It is similar to the `CMAKE_INSTALL_PREFIX` variable in CMAKE.

As an example, let me attach the full configuration command here.
```sh
$ ./configure --with-x=0 -with-pic --with-make-np=12 --with-mpi-compilers=1 --with-mpi-dir=/home/jliu/lib/mpich-3.3.2/ --with-scalar-type=real --with-precision=double --with-chaco=1 --download-chaco --with-hypre=1 --download-hypre --with-spai=1 --download-spai --with-sundials=1 --download-sundials --with-mumps=1 --download-mumps --with-scalapack=1 --download-scalapack --with-blacs=1 --download-blacs --with-spooles=1 --download-spooles --with-superlu_dist=1 --download-superlu_dist --with-superlu=1 --download-superlu --download-fblaslapack --download-parmetis --with-ml=1 --download-ml --with-eigen=1 --download-eigen --with-debugging=no COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native" --prefix=~/lib/petsc-3.8.4-opt --download-metis
```
