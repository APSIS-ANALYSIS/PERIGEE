# Configuration setup for Stampede2

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /work/01346/liujuy/stampede2/lib/VTK-7.1.1/lib/cmake/vtk-7.1)

set(PETSC_DIR /work/01346/liujuy/stampede2/lib/petsc-3.6.4-opt)

set(PETSC_ARCH .)

set(HDF5_ROOT /opt/apps/intel17/impi17_0/phdf5/1.8.16/x86_64 ) 

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
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ===================================================================
# Compiler options
# ===================================================================
set(CMAKE_C_COMPILER  /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicc)
set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-xCORE-AVX512 -O3 -xhost -Wall")
set(CMAKE_BUILD_TYPE RELEASE)
set( CMAKE_VERBOSE_MAKEFILE OFF )
set(CMAKE_CXX_STANDARD 11)

# END OF FILE
