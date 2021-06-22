# External Library Install Guide

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
[![made-with-bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)
![GitHub](https://img.shields.io/github/license/ju-liu/PERIGEE)

This is a guide for installing the VTK, PETSc, and related libraries, which are needed by the PERIGEE code. The default operating system is Linux Ubuntu. Installing these libraries on Mac OS may require minor changes in the following steps.

In general, we do an out-of-source compiling, except for PETSc. We will have all the libraries installed in a specified location (in this guide, $HOME/lib). This is because we typically do not have the access to /usr/local in clusters. Essentially, this requires you to specify the install prefix in the configuration stage for each library build.

Some of the libraries are not required for building the code. Typically, we require one to have VTK, PETSc, HDF5, and METIS as minimum requirement.

## Table of Contents

- [Create a lib folder](#Create-a-lib-folder)
- [Install Valgrind](#Install-Valgrind)
- [Install CMake](#Install-CMake)
- [Install MPICH](#Install-MPICH)
- [Install VTK](#Install-VTK)
- [Install PETSc](#Install-PETSc)
- [Install HDF5](#Install-HDF5)
- [Install METIS](#Install-METIS)
- [Install Gmsh](#Install-Gmsh)
- [Install ParaView](#Install-ParaView)

## Create a lib folder
It is recommended to have all the libraries in a single lib folder. Now we create an empty folder in the $HOME director.

```sh
$ cd $HOME
$ mkdir lib
```

## Install Valgrind
The Valgrind is a software for detecting memory leaking errors. It is **not mandatory** for compiling the PERIGEE code. For example, we do not install it on clusters. It is just useful when one do code development. Install Valgrind can be very straightforward in Linux.
```sh
$ sudo apt install valgrind
```
Then try
```sh
$ valgrind ls -l
```
to make sure it is working.

Of course, for experienced users, you can download the Valgrind source file and install it in the lib folder. I do not intend to cover that in this guide.


## Install CMake
Make sure you have cmake installed in your system. In case you do not have one, do the following.
```sh
$ sudo apt-get install cmake
```

If you are an advanced user and wants to build cmake youself, consult this page: https://cmake.org/install/

If you are using a cluster, cmake can be loaded by
```sh
$ module load cmake
```
## Install MPICH
First, download the source file, extract the tar bar, and rename the folder as a source folder:
```
$wget https://www.mpich.org/static/downloads/3.3rc1/mpich-3.3rc1.tar.gz
tar -zxvf mpich-3.3rc1.tar.gz
mv mpich-3.3rc1 mpich-3.3rc1-src
```
Note that you may want to go to https://www.mpich.org/static/downloads/ to select the version of the MPICH implementation. Versions 3.2 and 3.3 are conservative choices. Now you may enter the folder and do the following to install MPICH at the prescribed location.
```
./configure --prefix=$HOME/lib/mpich-3.3rc1 2>&1 | tee c.txt
make 2>&1 | tee m.txt
make install 2>&1 | tee mi.txt
```
You may want to add the bin folder to PATH:
```
vi ~/.bashrc
export PATH=mpich-install/bin:$PATH
```
so that the system will call the installed mpich binaries.


## Install VTK
The PERIGEE code is compatible with VTK-5, VTK-6, and VTK-7. In the following, we demonstrate steps for compiling VTK-7.1.1.

First, download the source file of VTK-7.1.1, and extract the tar ball. You do not have to download it to the lib folder because eventually we will remove the source folder afte we build the library. In this guide, we download it to the home directory.
```sh
$ cd $HOME
$ wget https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz
$ tar -zxvf VTK-7.1.1.tar.gz 
```
You will have a folder named VTK-7.1.1 now. To avoid confusion, we rename it as VTK-7.1.1-src to indicate this is a source file folder.

```sh
$ mv VTK-7.1.1 VTK-7.1.1-src
```
Second, create an out-of-source build folder, say, in your home directory. This folder will also be removed after building the library. Run the following commands.

```sh
$ cd $HOME
$ mkdir build_vtk
$ cd build_vtk
$ cmake $HOME/VTK-7.1.1-src -DCMAKE_INSTALL_PREFIX=$HOME/lib/VTK-7.1.1-shared -DBUILD_SHARED_LIBS=ON
$ make -j 6
$ make install
```
A folder named VTK-7.1.1-shared will be created in `$HOME/lib`. The name of this folder indicate that the library is a shared library, because we set `BUILD_SHARED_LIBS=ON` in the configuration stage. For an explanation of the differences between shared libraries and static libraries, please consult this [page](https://stackoverflow.com/questions/2649334/difference-between-static-and-shared-libraries). Now, if you go to `$HOME/lib/VTK-7.1.1-shared`, you should be able to see four sub-directories named as share, include, lib, and bin. Now, you may safely remove VTK-7.1.1.tar.gz, VTK-7.1.1-src, and build_vtk.

```sh
$ rm -rf VTK-7.1.1.tar.gz VTK-7.1.1-src build_vtk
```
For more detailed discussion of installing the VTK library, refer to the [guide for advanced users](install-advanced.md).

## Install PETSc
First, we enter the lib folder and download petsc source files in the folder.
```sh
$ cd $HOME/lib
$ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.3.tar.gz
$ tar -zxvf petsc-lite-3.11.3.tar.gz
$ ./configure
```
Follow the instructions from the configure command output. You may need to download necessary packages by appending argument to the configure command, and these arguments are system-dependent. For most systems, it is safe to do the following.
```sh
$ ./configure --download-mpich --download-fblaslapack
```
You may have more arguments for the configure command to control your installation. For more information, refer to the [advanced user guide](install-advanced.md).
If the configuration stage is completed successfully, there will be messages to guide you to complete the installation. As an example, on my own machine, I have to do the following steps.
```sh
$ make PETSC_DIR=/home/jliu/lib/petsc-3.11.3 PETSC_ARCH=arch-linux2-c-debug all
$ make PETSC_DIR=/home/jliu/lib/petsc-3.11.3 PETSC_ARCH=arch-linux2-c-debug check
```
On TACC machines, you do not have to install PETSc yourself. Run the following command, you will see the location of the PETSc library on TACC machines.
```sh
$ module load petsc
$ echo $TACC_PETSC_DIR
```

## Install HDF5
Run the following commands to install HDF5.
```sh
$ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.16/src/hdf5-1.8.16.tar.gz
$ tar -zxvf hdf5-1.8.16.tar.gz 
$ mv hdf5-1.8.16 hdf5-1.8.16-src
$ cd hdf5-1.8.16-src
$ ./configure --prefix=$HOME/lib/hdf5-1.8.16 --enable-cxx
$ make
$ make check
$ make install
$ make check-install
$ cd ..
$ rm -rf hdf5-1.8.16-src
```
In the lib folder, there should be a sub-directory hdf5-1.8.16, containing the HDF5 library. Again, on TACC machines, the HDF5 library is already installed. To locate it, run the following.
```sh
$ module load hdf5
$ echo $TACC_HDF5_DIR
```

## Install METIS
Run the following commands to install METIS
```sh
$ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-5.0.3.tar.gz
$ tar -zxvf metis-5.0.3.tar.gz
$ mv metis-5.0.3 metis-5.0.3-src
$ cd metis-5.0.3-src
$ make config prefix=$HOME/lib/metis-5.0.3
$ make
$ make install
$ cd ..
$ rm -rf metis-5.0.3-src
```

## Install Gmsh
The PERIGEE code does **not** depend on Gmsh. However, we frequently use Gmsh to generate input files for unstructured grid. Therefore, it is recommended that you install Gmsh on your machine. This step is **unnecessary** for compiling PERIGEE on clusters. Currently, the code is compatible with the mesh format from Gmsh-3. Run the following commands to install Gmsh. We do not need to link the PERIGEE code to Gmsh. For many cases, we generate a mesh using Gmsh as a binary executable.
```sh
$ wget http://gmsh.info/src/gmsh-3.0.6-source.tgz
$ tar -zxvf gmsh-3.0.6-source.tgz
$ cd gmsh-3.0.6-source
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/lib/gmsh-3.0.6
$ make
$ make install
$ cd ../..
$ rm -rf gmsh-3.0.6-source
```

## Install ParaView
If you really want to try, you may download the ParaView source file and do an out-of-source build, just like what we did for VTK. I just directly download the binary file from the [ParaView website](https://www.paraview.org/download/), and extract it in the `$HOME/lib` folder. Note, you do **not** have to install ParaView on a cluster.
