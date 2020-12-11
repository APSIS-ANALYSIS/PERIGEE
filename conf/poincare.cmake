# Configuration setup for machine poincare mac laptop

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /Users/juliu/lib/VTK-8.2.0-shared/lib/cmake/vtk-8.2)

set(PETSC_DIR /Users/juliu/lib/petsc-3.11.3-debug)

set(PETSC_ARCH .)

set(HDF5_ROOT /Users/juliu/lib/hdf5-1.8.16)

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
set(CMAKE_C_COMPILER  /Users/juliu/lib/mpich-3.3.2/bin/mpicc)
set(CMAKE_CXX_COMPILER  /Users/juliu/lib/mpich-3.3.2/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O0 -Wall")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
