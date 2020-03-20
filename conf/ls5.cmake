# Configuration setup for Lonestar5 at TACC

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /work/01346/liujuy/lonestar/lib/VTK-6.3.0/lib/cmake/vtk-6.3)
#set(PETSC_DIR /work/01346/liujuy/lonestar/lib/petsc-3.8.4-opt-largemem)
#set(PETSC_ARCH .)

set(PETSC_DIR /opt/apps/intel18/cray_mpich_7_7/petsc/3.11)
set(PETSC_ARCH haswell)

#set(HDF5_DIR /opt/apps/intel18/hdf5/1.8.16/x86_64)
set(HDF5_ROOT ${PETSC_DIR}/${PETSC_ARCH})

#SET(SLEPC_DIR /work/01346/liujuy/lonestar/lib/slepc-3.10.1-opt )

#SET(BOOST_DIR /opt/apps/intel18/boost/1.64)

#SET(BESSEL_DIR /work/01346/liujuy/lonestar/lib/complex-bessel)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)
#find_package(SLEPc)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})
include_directories(${HDF5_INCLUDE_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ===================================================================
# Compiler options
# ===================================================================
set(CMAKE_C_COMPILER  /opt/apps/intel18/cray_mpich/7.7.3/bin/mpicc)
set(CMAKE_CXX_COMPILER  /opt/apps/intel18/cray_mpich/7.7.3/bin/mpicxx)
#set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/mpicc)
#set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O3 -xhost -Wall")
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
