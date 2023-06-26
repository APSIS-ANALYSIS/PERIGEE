# Configuration setup for machine kolmogorov mac laptop

# ========================================================
# Specify the library locations
# ========================================================
set(HOME_DIR /Users/juliu)

set(VTK_DIR ${HOME_DIR}/lib/VTK-9.2.6-SHARED/lib/cmake/vtk-9.2)

set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.18.6-debug)
set(SLEPC_DIR ${HOME_DIR}/lib/slepc-3.18.3-debug)

set(PETSC_ARCH .)

set(HDF5_ROOT ${HOME_DIR}/lib/hdf5-1.12.2)

# ==========================WOMERSLY CHANGE==========================
set(BESSEL_DIR /Users/juliu/lib/complex-bessel)
# ==========================WOMERSLY CHANGE==========================

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK COMPONENTS CommonCore CommonSystem
  CommonMisc CommonMath IOCore IOLegacy IOXML REQUIRED)
find_package(HDF5 REQUIRED)
find_package(SLEPc)

set(ENV{PKG_CONFIG_PATH} ${PETSC_DIR}/lib/pkgconfig)

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${SLEPC_INC})
include_directories(${PETSC_INCLUDE_DIRS})

# ==========================WOMERSLY CHANGE==========================
include_directories(${BESSEL_DIR}/include)
include_directories(${BESSEL_DIR}/include/complex_bessel_bits)
# ==========================WOMERSLY CHANGE==========================


link_directories(${PETSC_STATIC_LIBRARY_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SLEPC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_STATIC_LIBRARIES})

# ==========================WOMERSLY CHANGE==========================
LINK_DIRECTORIES( ${BESSEL_DIR}/lib )
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} complex_bessel)
# ==========================WOMERSLY CHANGE==========================

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER ${PETSC_DIR}/bin/mpicc)
set(CMAKE_CXX_COMPILER ${PETSC_DIR}/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
