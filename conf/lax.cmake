# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

set(VTK_DIR /home/jliu/lib/VTK-8.1.1/lib/cmake/vtk-8.1)

if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-opt)
  set(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-opt)
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-debug)
  set(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-debug)
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(PETSC_ARCH .)

set(METIS_DIR /home/jliu/lib/metis-5.0.3)

SET(HDF5_ROOT /home/jliu/lib/hdf5-1.8.20)

#include_directories(/home/jliu/lib/eigen-3.3.4/include/eigen3)

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

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicc)
set(CMAKE_CXX_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicxx)
set(CMAKE_CXX_STANDARD 11)
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O3 -Wall")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

message(STATUS ${EXTRA_LINK_LIBS})

# This an example for buiding code with prof for profiling
#SET(CMAKE_CXX_FLAGS "-pg -Wall")
