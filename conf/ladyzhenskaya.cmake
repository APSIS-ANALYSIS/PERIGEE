# Configuration setup for machine ladyzhenskaya Linux desktop

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /home/ju/lib/VTK-8.2.0-debug/lib/cmake/vtk-8.2)

if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/ju/lib/petsc-3.14.2-opt)
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/ju/lib/petsc-3.14.2-debug)
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(PETSC_ARCH .)

set(HDF5_ROOT /home/ju/lib/hdf5-1.8.16)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy vtkIOXML REQUIRED)
find_package(HDF5 REQUIRED)
find_package(PETSc REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PETSC_INC})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  /home/ju/lib/mpich-3.3.2/bin/mpicc)
set(CMAKE_CXX_COMPILER  /home/ju/lib/mpich-3.3.2/bin/mpicxx)
set(CMAKE_CXX_STANDARD 11)

if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O3 -Wall")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
