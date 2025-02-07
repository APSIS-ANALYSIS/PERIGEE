## Install PERIGEE on Windows

To install PERIGEE on Windows, we must first install Windows Subsystem for Linux (WSL). Then, we may perform procedures similar to those of UNIX-like systems to complete the install.

### Install WSL

For details, refer to the [official guidance](https://learn.microsoft.com/zh-cn/windows/wsl/). According to the official guidance, WSL can be installed simply by running the following in "Command Prompt" as administrator and restarting on the Windows with OS build later than [19041](https://en.wikipedia.org/wiki/Windows_10_version_history).
```sh
wsl --install -d Ubuntu-20.04
```

Here, I post my steps to install WSL.

Open "Control Panel > Programs > Turn Windows features on or off". Select the checkbox of "Windows Subsystem for Linux", click "OK", and then restart your Windows.

Run "Command Prompt" as administrator. Run the following to check the valid distributions that can be installed.
```sh
wsl -l -o
```

Then, run the following to install a distribution. Here, "Ubuntu 20.04 LTS" is recommended.
```sh
wsl --install -d Ubuntu-20.04
```

Then, Open the application "Ubuntu 20.04.6 LTS". It will ask you to enter a UNIX name and set a password. After that, you have the WSL successfully installed. If it throws "The attempted operation is not supported for the type of object referenced." Open "Command Prompt" as administrator. Run the following and restart your WSL, not the Windows.
```sh
net winsock reset
```
You can search the following in "This PC" to access WSL folder.
```sh
\\wsl$
```
Additionally, you can access your Windows disk by running the following in WSL.
```sh
cd /mnt/c
```
You can run the following commands to set up your WSL.
```sh
sudo apt update
sudo apt upgrade
sudo apt autoremove
sudo apt install vim
sudo apt install git-all
sudo apt install build-essential
sudo apt install gfortran
sudo apt install python2
sudo apt install python3
sudo apt install mesa-utils
sudo apt install mesa-common-dev
sudo apt install libgl1-mesa-dev
sudo apt install libxt-dev
sudo apt install cmake
```

### Install External Libs
Create a lib directory by running:
```sh
cd $HOME
mkdir lib
```
#### 1. Install MPICH
Go to the lib directory.
```sh
cd $HOME/lib
```
Download the source, extract the tar bar, and rename the directory as a source directory.
```sh
wget https://www.mpich.org/static/downloads/3.3rc1/mpich-3.3rc1.tar.gz
tar -zxvf mpich-3.3rc1.tar.gz
mv mpich-3.3rc1 mpich-3.3rc1-src
```
Enter the folder and install MPICH.
```sh
cd mpich-3.3rc1-src
./configure --prefix=$HOME/lib/mpich-3.3rc1 2>&1 | tee c.txt
```
When the configuration is complete, run the following.
```sh
make 2>&1 | tee m.txt
make install 2>&1 | tee mi.txt
```
After that, you can remove the useless files and directories.
```sh
cd $HOME/lib
rm -r mpich-3.3rc1-src mpich-3.3rc1.tar.gz
```
Then, add the bin directory to the system variable PATH.
```sh
vi ~/.bashrc
```
Add the following to .bashrc
```sh
export PATH=$HOME/lib/mpich-3.3rc1/bin:$PATH
```
Exit ``.bashrc`` and source it, the system will call the installed MPICH.
```sh
source ~/.bashrc
```

#### 2. Install VTK
Go to the lib directory, download VTK source, extract the tar bar, and rename the directory as a source directory.
```sh
cd $HOME/lib
wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar -zxvf VTK-8.2.0.tar.gz
mv VTK-8.2.0 VTK-8.2.0-src
```
Create a build directory and run the cmake to compile VTK.
```sh
cd $HOME/lib
mkdir build_vtk
cd build_vtk
cmake ../VTK-8.2.0-src/ -DCMAKE_INSTALL_PREFIX=$HOME/lib/VTK-8.2.0-shared -DBUILD_SHARED_LIBS=ON  -DCMAKE_BUILD_TYPE=Release
make -j 4
make install
```
Remove the useless files when VTK is successfully installed.
```sh
cd $HOME/lib
rm -rf build_vtk VTK-8.2.0-src VTK-8.2.0.tar.gz
```
After that, you need to add the VTK path to your ``.bashrc``.
```sh
export LD_LIBRARY_PATH=/home/yxh/lib/VTK-8.2.0-shared/lib:$LD_LIBRARY_PATH
```
#### 3. Install PETSc
Go to the lib directory, download PETSc source, extract the tar bar, and rename the directory as a source directory.
```sh
cd $HOME/lib
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.18.6.tar.gz
tar -zxvf petsc-3.18.6.tar.gz
mv petsc-3.18.6 petsc-3.18.6-src
cd petsc-3.18.6-src
```
Then, we need to run the configuration. It is recommended to write a script.
```sh
vi conf.sh
```
Copy all the following to the script file, ``conf.sh``.
```sh
./configure --with-x=0 \
  -with-pic \
  --with-make-np=4 \
  --with-mpi-compilers=1 \
  --with-mpi-dir=${HOME}/lib/mpich-3.3rc1/ \
  --with-scalar-type=real \
  --with-precision=double \
  --with-hypre=1 \
  --download-hypre \
  --with-mumps=1 \
  --download-mumps \
  --with-scalapack=1 \
  --download-scalapack \
  --with-blacs=1 \
  --download-blacs \
  --download-fblaslapack \
  --download-metis \
  --download-hdf5 \
  --with-debugging=no \
  COPTFLAGS="-O3 -march=native -mtune=native" \
  CXXOPTFLAGS="-O3 -march=native -mtune=native" \
  FOPTFLAGS="-O3 -march=native -mtune=native" \
  --prefix=${HOME}/lib/petsc-3.18.6-opt
```
Save and exit the script file, run it by
```sh
sh conf.sh
```
After the configuration, run the following as the system asked
```sh
make PETSC_DIR=/home/yxh/lib/petsc-3.18.6-src PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=/home/yxh/lib/petsc-3.18.6-src PETSC_ARCH=arch-linux-c-opt install
make PETSC_DIR=/home/yxh/lib/petsc-3.18.6-opt PETSC_ARCH="" check
```
Then, remove the useless directory and file.
```sh
cd $HOME/lib
rm -rf petsc-3.18.6-src petsc-3.18.6.tar.gz
```

### Install PERIGEE
First, download the source file by git.
```sh
cd $HOME
mkdir code
cd code
git clone https://github.com/APSIS-ANALYSIS/PERIGEE.git
```
Then, edit the configuration file.
```sh
cd PERIGEE/conf
vi system_lib_loading.cmake
```
Copy the following to your file, and carefully set your home directory.
```sh
# Configuration setup for machine kolmogorov mac laptop
# ========================================================
# Specify the library locations
# ========================================================
set(HOME_DIR /home/yxh) 
set(VTK_DIR ${HOME_DIR}/lib/VTK-8.2.0-shared/lib/cmake/vtk-8.2)
set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.18.6-opt)
set(PETSC_ARCH .)
set(METIS_DIR ${PETSC_DIR})
set(HDF5_ROOT ${PETSC_DIR})
# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(VTK REQUIRED)
find_package(HDF5 REQUIRED)
find_package(PETSc REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
message(STATUS "Use METIS in PETSC: " ${PETSC_METIS_LIB})
message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})
# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER ${HOME_DIR}/lib/mpich-3.3rc1/bin/mpicc)
set(CMAKE_CXX_COMPILER ${HOME_DIR}/lib/mpich-3.3rc1/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE OFF)
# EOF
```
Save and exit the file. Next, make a build directory.
```sh
cd $HOME
mkdir build
cd build
```
Run the following to build the code.
```sh
cmake $HOME/code/PERIGEE/examples/ruc_fsi/ -DCMAKE_BUILD_TYPE=Release
make -j 4
```
