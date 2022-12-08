# PERIGEE Configuration Guide

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
[![made-with-bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)
![GitHub](https://img.shields.io/github/license/ju-liu/PERIGEE)

In this guide, we assume the external libraries have been correctly installed, following the steps in the [installation guide](install-external-libs.md). We will essentially tell the PERIGEE code where are those libraries.

## Table of Contents
- [Name the machine](#Name-the-machine)
- [Modify the CMake file](#Modify-the-Cmake-file)
- [Create another CMake file for a cluster (optional)](#Create-another-CMake-file-for-a-cluster-(optional))
- [Modify the system_lib_loading.cmake file](#Modify-the-system_lib_loading.cmake-file)

## Name your machine
Before you start, I recommend giving a meaningful name for your machine. What I typically do is that I give an environmental variable for my machine to name it. My first office desktop is called Bacon, and I have been using this name since then. To give your machine a name, go back to the home directory and open the .bashrc file, add `export MACHINE_NAME=bacon` to this file. Save it and close the terminal. Reopen the terminal and type
```sh
$ echo $MACHINE_NAME
```
You shall see bacon on screen. With a manchine name, the configuration may detect the machine, and load the corresponding cmake file. This will make the configuration across platforms easier. For TACC machines, this environmental variable is already defined by them. For example, in Stampede2, the value of `$MACHINE_NME` is `sp2`.

## Modify the CMake file
There is a folder named conf in the PERIGEE code. It stores the CMake configuration files for compiling the code. You may find a file called `machine_name.cmake`. Now, let us first make a copy, with a different meaningful name.
```sh
$ cp machine_name.cmake linux_bacon.cmake
```
Now we open the linux_bacon.cmake and there will be a few variables to modify.

1. Set `VTK_DIR` to the VTK library location (e.g. `/home/jliu/lib/VTK-7.1.1-shared`).
2. Set `PETSC_DIR` to the PETSc library location (e.g. `/home/jliu/lib/petsc-3.11.3`).
3. Set `PETSC_ARCH` to the value used in PETSc installation (e.g. `arch-linux2-c-debug`).
4. Set `METIS_DIR` to the METIS library location (e.g. `/home/jliu/lib/metis-5.0.3`).
5. Set `HDF5_DIR` to the HDF5 library location (e.g. `/home/jliu/lib/hdf5-1.8.16`).
6. Set `CMAKE_C_COMPILER` to `$PETSC_DIR/$PETSC_ARCH/bin/mpicc`
7. Set `CMAKE_CXX_COMPILER` to `$PETSC_DIR/$PETSC_ARCH/bin/mpicxx`

After the editing, save and close linux_bacon.cmake. And you have your own configuration file setted up.

## Create another CMake file for a cluster (optional)
One is highly likely to edit, test, and debug the code on a local desktop, and the production code will be run on a cluster. So one may have to create another CMake configuration file for the cluster. There is a file named stampede2_juliu.cmake in the conf folder. It is quite similar to the linux_bacon.cmake file with different settings of the CMake variables. You are recommended to do the following to create your own CMake file first.
```sh
$ cp machine_name.cmake stampede2_yourname.cmake
```
And modify the file. Of course, you may need to consult [installation guide](install-external-libs.md) to build necessary libraries. On Stampede2, PETSc and HDF5 are readily available. You need to have VTK, METIS, and (optionally) Gmsh installed.

## Modify the system_lib_loading.cmake file
Now, we may open system_lib_loading.cmake in conf folder. The file contains a two-level if-elseif statement to decide which .cmake configuration file to load. On the first level, the if-statement detects if the machine is a Mac OS or Linux system. I have a Mac Laptop, and you may see that my laptop is called poincare, and its cmake file is `poincare_PETSc_VTK.cmake`. On the second level, there is an if-elseif statement to check the system's environmental variable `$ENV{MACHINE_NAME}`. This is the variable you edited in the .bashrc file in the beginning of this guide. You may follow the code style add an elseif statement to load your cmake file.
