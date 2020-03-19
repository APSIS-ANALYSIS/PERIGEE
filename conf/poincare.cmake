# Configuration setup for machine poincare mac laptop

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /Users/juliu/lib/VTK-7.1.1-debug/lib/cmake/vtk-7.1)

set(PETSC_DIR /Users/juliu/lib/petsc-3.11.3-debug)

set(PETSC_ARCH .)

set(HDF5_ROOT ${PETSC_DIR})

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})

include_directories(${HDF5_INCLUDE_DIRS})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  /Users/juliu/lib/mpich-3.2.1/bin/mpicc)
set(CMAKE_CXX_COMPILER  /Users/juliu/lib/mpich-3.2.1/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O0 -Wall -Winconsistent-missing-override")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
