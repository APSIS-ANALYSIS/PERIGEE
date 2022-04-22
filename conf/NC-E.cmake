# Configuration setup for NC-E ningxia

set(HOME_DIR /public1/home/scfa1548/)

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR ${HOME_DIR}/lib/VTK-7.1.1-opt/lib/cmake/vtk-7.1)

# gcc
# set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.15.5-mpich3.2-gcc5.4.0-opt)

# icc
set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.13.6-opt)

set(PETSC_ARCH .)

set(HDF5_ROOT ${HOME_DIR}/lib/hdf5-1.12.1)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})
include_directories(${HDF5_INCLUDE_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

if(PETSC_METIS)
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
  message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
else(PETSC_METIS)
  find_package(METIS)
  include_directories(${METIS_INCLUDE_DIRS})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LIBRARIES})
endif(PETSC_METIS)

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ===================================================================
# Compiler options
# ===================================================================
set(CMAKE_C_COMPILER /public1/soft/intel/2018.4/compilers_and_libraries/linux/mpi/intel64/bin/mpicc)
set(CMAKE_CXX_COMPILER /public1/soft/intel/2018.4/compilers_and_libraries/linux/mpi/intel64/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF

