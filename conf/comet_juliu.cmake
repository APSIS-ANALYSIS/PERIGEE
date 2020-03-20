# Configuration setup for machine comet

# Library locations
set(VTK_DIR /oasis/scratch/comet/liujuy/temp_project/lib/VTK-7.1.1/lib/cmake/vtk-7.1)

set(PETSC_DIR /oasis/scratch/comet/liujuy/temp_project/lib/petsc-3.7.7-opt )

set(PETSC_ARCH . )

set(HDF5_ROOT /opt/hdf5/1.8.21/intel/mvapich2_ib )

# Setup the libraries
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})
include_directories(${HDF5_INCLUDE_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# 5. Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  ${PETSC_DIR}/bin/mpicc)
set(CMAKE_CXX_COMPILER ${PETSC_DIR}/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_CXX_STANDARD 11)

# EOF
