# Configuration setup for machine Lax

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /home/jliu/lib/VTK-8.1.1/lib/cmake/vtk-8.1)

if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-opt)
  set(PETSC_ARCH .)
  set(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-opt)
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-debug)
  set(PETSC_ARCH .)
  set(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-debug)
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(METIS_DIR /home/jliu/lib/metis-5.0.3)

set(HDF5_ROOT /home/jliu/lib/hdf5-1.8.20)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)
find_package(SLEPc)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})

if(PETSC_METIS)
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
  message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
else(PETSC_METIS)
  find_package(METIS)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LIBRARIES})
endif(PETSC_METIS)

include_directories(${HDF5_INCLUDE_DIRS})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

if(SLEPC_FOUND)
  include_directories(${SLEPC_INC})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SLEPC_LIB})
endif(SLEPC_FOUND)

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_CXX_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicxx)
set(CMAKE_C_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicc)
set(CMAKE_CXX_STANDARD 11)
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O3 -Wall")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set( CMAKE_VERBOSE_MAKEFILE OFF )

# This an example for buiding code with prof for profiling
#SET(CMAKE_CXX_FLAGS "-pg -Wall")

# EOF
