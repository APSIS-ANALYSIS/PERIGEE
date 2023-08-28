# External Library Install Guide

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
[![made-with-bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)
![GitHub](https://img.shields.io/github/license/ju-liu/PERIGEE)

This is a guide for installing the libraries that PERIGEE is dependent on. The default operating system is Linux Ubuntu. Installing these libraries on Mac OS requires minor changes in the following steps.

You will have all the libraries installed in a specified location (in this guide, $HOME/lib). This is because one typically does not have an access to /usr/local in remote clusters. Roughly speaking, this requires you to specify the install prefix in the configuration stage for the build of each library.

The four libraries VTK, PETSc, HDF5, and METIS are required for the build of PERIGEE. MPICH is needed for parallelization; Gmsh is often needed for mesh generation; SLEPc is needed for calculating eigenvalues; ParaView is used for visualization solution outputs.

On most clusters, there are already certain libraries installed by the admin. You need to use [Environment Modules](https://modules.readthedocs.io/en/latest/) to manage the installed libraries if you want to use them. As a short note, you can use `module list` to see the currently loaded modules and `module avail` to see the available (or, installed) modules. You may load and unload one by `module load [name]` and `module unload [name]`. You may use `module purge` to unload all modules. If you want to see the detailed info about the module, in particular the library path, you may use `module show [name]` to list related info. As an example, in my `.bashrc` file of TaiYi machine, I have used the following.
```
module purge
module load intel/2018.4
module load mpi/intel/2018.4
```
The above indicates that all modules are unloaded first and the intel/2018.4 compiler is loaded, and then the corresponding MPI is loaded.

## Table of Contents
- [Setup the machine](#Setup-the-machine)
- [Create a lib folder](#Create-a-lib-folder)
- [Install Valgrind](#Install-Valgrind)
- [Install CMake](#Install-CMake)
- [Install MPICH](#Install-MPICH)
- [Install VTK](#Install-VTK)
- [Install PETSc](#Install-PETSc)
- [Install HDF5](#Install-HDF5)
- [Install METIS](#Install-METIS)
- [Install Gmsh](#Install-Gmsh)
- [Install SLEPc](#Install-SLEPc)
- [Install ParaView](#Install-ParaView)

## Setup the machine
Assuming that you just installed a Ubuntu system. You may want to do the following to setup the system with essiential components.
```sh
sudo apt update
sudo apt upgrade
sudo apt autoremove
sudo apt install vim
sudo apt install git-all
sudo apt install build-essential
sudo apt install texlive-latex-extra
sudo apt install texlive-publishers
sudo apt install texlive-science
sudo apt install gfortran
sudo apt install python2
sudo apt install python3
sudo apt install mesa-utils
sudo apt install mesa-common-dev
sudo apt install libgl1-mesa-dev
sudo apt install libxt-dev
sudo apt install cmake
sudo apt install valgrind
```
Go to Ubuntu Software and install Texmaker, HDFView, and JabRef.

For Acrobat Reader, you may need to run the following.
```sh
sudo apt install gdebi-core libxml2:i386 libcanberra-gtk-module:i386 gtk2-engines-murrine:i386 libatk-adaptor:i386
wget ftp://ftp.adobe.com/pub/adobe/reader/unix/9.x/9.5.5/enu/AdbeRdr9.5.5-1_i386linux_enu.deb
sudo gdebi AdbeRdr9.5.5-1_i386linux_enu.deb
```
Then clean the .deb file from your disk.

For machines with NVIDA graphic processors, you may need to install a driver. Run `ubuntu-drivers devices` to see the type of devices. Run `sudo ubuntu-drivers autoinstall` to install the recommended driver. Or, you may install specific driver by `sudo apt install nvidia-340`. You may refer to this [website](https://linuxconfig.org/benchmark-your-graphics-card-on-linux) to test your graphic card. For more detailed guidance, you may refer to this [guide](https://zhuanlan.zhihu.com/p/59618999).

## Create a lib folder
It is recommended to have all the libraries in a single lib folder. We typically create an empty folder in the $HOME directory.

```sh
$ cd $HOME
$ mkdir lib
```

## Install Valgrind
The Valgrind is a software for detecting memory leaking errors. It is **not mandatory** for compiling the PERIGEE code. For example, we do not install it on clusters. It is just useful in code development. Installing Valgrind can be very straightforward in Linux.
```sh
$ sudo apt install valgrind
```
Then try
```sh
$ valgrind ls -l
```
to make sure it is working.

Of course, for experienced users, you can download the Valgrind source file and install it in the lib folder.


## Install CMake
Make sure you have cmake installed in your system. In case you do not have one, you will have to install it youself. There are several ways of intalling CMake. If you have the admin privilege, do the following.
```sh
$ sudo apt-get install cmake
```
Linux will automatically install CMake for you. If you are an advanced user and wants to build cmake youself, consult this page: https://cmake.org/install/ You may encounter issues like the lack of openssl in the system, which can be conveniently resolved by, e.g., `sudo apt install libssl-dev`. Of course to tell the system to load the correct CMake, you may also want to modify PATH in the .bashrc or .bash_profile file.

If you are using a cluster, cmake can be loaded by
```sh
$ module load cmake
```
The easiest way of getting CMake is directly downloading a compiled CMake from the [official site](https://cmake.org/download/). When you download the compiled CMake, make sure you get the correct platform. Typically, it should be Linux x86-64.

## Install MPICH
First, download the source file, extract the tar bar, and rename the folder as a source folder:
```
$ wget https://www.mpich.org/static/downloads/3.3rc1/mpich-3.3rc1.tar.gz
$ tar -zxvf mpich-3.3rc1.tar.gz
$ mv mpich-3.3rc1 mpich-3.3rc1-src
```
Remark: you may want to go to https://www.mpich.org/static/downloads/ to select the version of the MPICH implementation. Versions 3.2 and 3.3 are conservative choices. Now you may enter the folder and do the following to install MPICH at the prescribed location.
```
./configure --prefix=$HOME/lib/mpich-3.3rc1 2>&1 | tee c.txt
```
At this stage, you may encounter a warning saying you do not have Fortran compilers. You may want to install gfortran by 
```
sudo apt-get install gfortran
```
and rerun the configure command. Or you may encounter an issue saying `gfortran will not compile files that call the same routine with arguments of different types`. You need to add the following to your bashrc or bash_profile file.
```
export FFLAGS="-w -fallow-argument-mismatch -O2"
```
Once you see that the system tells the configuration is complete, run the following.
```
make 2>&1 | tee m.txt
make install 2>&1 | tee mi.txt
```
You may want to add the bin directory to the system variable PATH. First, open the .bashrc file by vi.
```
vi ~/.bashrc
```
Add the following statement to .bashrc
```
export PATH=$HOME/lib/mpich-3.3rc1/bin:$PATH
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
```
At this stage, you may encounter an issue due to the lack of the OpenGL library. There are two ways to address this issue. Suppose the machine is of your own and you have the admin privilege. You may choose to install OpenGL. Otherwise, if the machine is a cluster and you are a regular user, you may install VTK without OpenGL. Let us give the first approach here. Run `sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev`. The second approach is documented in [guide for advanced users](install-advanced.md).
```
$ make -j 6
$ make install
```
A folder named VTK-7.1.1-shared will be created in `$HOME/lib`. The name of this folder indicate that the library is a shared library, because we set `BUILD_SHARED_LIBS=ON` in the configuration stage. For an explanation of the differences between shared libraries and static libraries, please consult this [page](https://stackoverflow.com/questions/2649334/difference-between-static-and-shared-libraries). Now, if you go to `$HOME/lib/VTK-7.1.1-shared`, you should be able to see four sub-directories named as share, include, lib, and bin. With the library installed, you may safely remove VTK-7.1.1.tar.gz, VTK-7.1.1-src, and build_vtk.

```sh
$ rm -rf VTK-7.1.1.tar.gz VTK-7.1.1-src build_vtk
```
For more detailed discussion of installing the VTK library, refer to the [guide for advanced users](install-advanced.md). We recommend using VTK-7.1.1 or VTK-8.2.0.

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

## Install HDF5
Run the following commands to install HDF5.
```sh
$ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.16/src/hdf5-1.8.16.tar.gz
$ tar -zxvf hdf5-1.8.16.tar.gz 
$ mv hdf5-1.8.16 hdf5-1.8.16-src
$ cd hdf5-1.8.16-src
$ ./configure --prefix=$HOME/lib/hdf5-1.8.16 --enable-production
```
So far, you have downloaded the source file and configured the HDF5 source. In the configuration stage, you specified the install location for the HDF5 library and the compile uses production mode rather than debug mode. You can always see all configuration options by running ```./configure --help```. At the end of configure output, there will be a summary of the configuration setup. Make sure that (1) the install location is correct; (2) the build type is production (or debug); (3) optionally, parallel hdf5 is supported.

```sh
$ make
$ make check
$ make install
$ make check-install
$ cd ..
$ rm -rf hdf5-1.8.16-src
```
The commands `make check` and `make check-install` are optional but recommended. In the lib folder, there should be a sub-directory hdf5-1.8.16, containing the HDF5 library. As a note, we recommend using hdf5-1.8.16 or hdf5-1.12.2.

## Install METIS
Run the following commands to install METIS.
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
Remarks: (1) You should read BUILD.txt to have an idea of all relevant configure options. (2) If you need 64-bit integer support, you need to modify the `metis.h` file in the include folder before compiling.

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

## Install SLEPc
The PERIGEE code in general does **not** depend on the SLEPc library. Yet, oen may occasionally need it to calculate eigenvalues. You may want to consult their [webpage](slepc.upv.es) for more information. Run the following to build the SLEPc library.
```sh
$ wget https://slepc.upv.es/download/distrib/slepc-3.11.3.tar.gz
$ tar -zxvf slepc-3.11.3.tar.gz
$ mv slepc-3.11.3 slepc-3.11.3-src
$ cd slepc-3.11.3-src
$ export PETSC_DIR=/home/juliu/lib/petsc-3.11.3-debug
$ export SLEPC_DIR=/home/juliu/lib/slepc-3.11.3-src
$ ./configure --prefix=/home/juliu/lib/slepc-3.11.3-debug
$ make SLEPC_DIR=/home/juliu/lib/slepc-3.11.3 PETSC_DIR=/home/juliu/lib/petsc-3.11.3-debug
$ make SLEPC_DIR=/home/juliu/lib/slepc-3.11.3 PETSC_DIR=/home/juliu/lib/petsc-3.11.3-debug install
$ make SLEPC_DIR=/home/juliu/lib/slepc-3.11.3-debug PETSC_DIR=/home/juliu/lib/petsc-3.11.3-debug PETSC_ARCH="" check
```
Notice that one needs to make the SLPEc's major version compatible with that of PETSc. Also, one may build an optimized library by setting the `PETSC_DIR` to the optimized PETSc build's path.

## Install ParaView
If you really want to try, you may download the ParaView source file and do an out-of-source build, just like what we did for VTK. I just directly download the binary file from the [ParaView website](https://www.paraview.org/download/), and extract it in the `$HOME/lib` folder. Note, you do **not** have to install ParaView on a cluster.
