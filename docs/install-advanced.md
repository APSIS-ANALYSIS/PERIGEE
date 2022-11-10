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
On certain clusters, OpenGL is not installed. So we have to turnoff the rendering. A minimum required install for VTK is as follows.
```sh
$ cmake ../VTK-7.1.1 -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_INSTALL_PREFIX:PATH=/home/juliu/lib/VTK-7.1.1-OPT -DVTK_Group_StandAlone:BOOL=OFF -DVTK_Group_Rendering:BOOL=OFF -DModule_vtkCommonMath:BOOL=ON -DModule_vtkCommonMisc:BOOL=ON -DModule_vtkCommonCore:BOOL=ON -DModule_vtkCommonSystem:BOOL=ON -DModule_vtkIOCore:BOOL=ON -DModule_vtkIOLegacy:BOOL=ON -DModule_vtkIOXML:BOOL=ON
$ make
$ make install
```
In the above, we only install the modules we need.

## Advanced Guide for PETSc Installation
The PETSc package has an official installation guide [page](https://www.mcs.anl.gov/petsc/documentation/installation.html). The package installation is controlled through the `configure` command.
* `--with-mpi-dir=/home/jliu/lib/mpich-3.2` The `--with-mpi-dir` tells the PETSc that there exists a MPI library installed in the computer, and the PETSc will not have to download and compile a MPICH during installation. My own mpich is installed in `/home/jliu/lib/mpich-3.2`.

* `--with-hypre=1 --download-hypre` This flag tells the PETSc installer to install the [Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) package as well. Several preconditioners in PETSc requires external libraries. We recommend one install at least Hypre and MUMPS in PETSc. 

Also, notice that sometimes, PETSc configure command may not have the correct download address for HYPRE. You may have to go to https://github.com/hypre-space/hypre to download the HYPRE source with the correct version manually. And then, add `--download-hypre=some address/hypre-source.tar.gz` in the configure option.

* `--with-hdf5=1 --download-hdf5 --with-metis=1 --download-metis` This tells PETSc to install HDF5 and METIS. With these two installed within PETSc, you will not need to install them separately youself. In the configuration file for PERIGEE, provide the correct path to link to them.

* `--with-debugging=yes` The `--with-debugging` flag is similar to the `CMAKE_BUILD_TYPE` in CMake. It tells the installer if you want to have a debug mode or an optimized mode for the library. Sometimes, the following will also be needed, `COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native"`.

* `--prefix=/home/jliu/lib/petsc-3.11.3-debug` The `--prefix` flag tells the installer where you want to install the PETSc package. It is similar to the `CMAKE_INSTALL_PREFIX` variable in CMAKE.

* Refer to [this](https://www.mcs.anl.gov/petsc/documentation/installation.html) for the installation details.

As an example, let me attach the full configuration command for debug build here.
```sh
 ./configure -with-pic --with-make-np=12 --with-mpi-compilers=1 --with-mpi-dir=/home/juliu/lib/mpich-3.2/ --with-scalar-type=real --with-precision=double --with-chaco=1 --download-chaco --with-hypre=1 --download-hypre=/home/juliu/lib/hypre-2.11.1.tar.gz --with-mumps=1 --download-mumps --with-scalapack=1 --download-scalapack --download-fblaslapack --download-metis --download-parmetis --with-yaml=1 --download-yaml --with-debugging=yes --prefix=/home/juliu/lib/petsc-xx.yy.zz-debug
```
Let me attach the full configuration command for optimized build here.
```sh
$ ./configure --with-x=0 -with-pic --with-make-np=12 --with-mpi-compilers=1 --with-mpi-dir=${HOME}/lib/mpich-3.3.2/ --with-scalar-type=real --with-precision=double --with-chaco=1 --download-chaco --with-hypre=1 --download-hypre --with-spai=1 --download-spai --with-sundials=1 --download-sundials --with-mumps=1 --download-mumps --with-scalapack=1 --download-scalapack --with-blacs=1 --download-blacs --with-spooles=1 --download-spooles --with-superlu_dist=1 --download-superlu_dist --with-superlu=1 --download-superlu --download-fblaslapack --download-metis --download-parmetis --with-ml=1 --download-ml --with-eigen=1 --download-eigen --with-yaml=1 --download-yaml --with-debugging=no COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native" --prefix=${HOME}/lib/petsc-xx.yy.zz-opt
```
In the above, you should specify the correct version number in `petsc-xx.yy.zz-opt`.

My configuration on TaiYi is as follows.
```sh
./configure --with-mpi-dir=/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/ --with-blaslapack-dir=/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mkl --with-debugging=no --prefix=/work/mae-liuj/lib/petsc-3.16.6-opt --download-hypre=/work/mae-liuj/petsc-3.16.6-extlibs/hypre-2.23.0.tar.gz --download-mumps=/work/mae-liuj/petsc-3.16.6-extlibs/petsc-pkg-mumps-6d1470374d32.tar.gz --download-metis=/work/mae-liuj/petsc-3.16.6-extlibs/petsc-pkg-metis-c8d2dc1e751e.tar.gz COPTFLAGS="-O3 -xHOST" CXXOPTFLAGS="-O3 -xHOST" FOPTFLAGS="-O3 -xHOST" --with-scalapack-include=/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mkl/include --with-scalapack-lib="-L/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64/ -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64"
```
