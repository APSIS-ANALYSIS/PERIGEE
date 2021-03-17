# Configuration setup for the Linux machine of Ingrid's

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /home/ingridxlan/lib/VTK-7.1.1-shared/lib/cmake/vtk-7.1)

set(PETSC_DIR /home/ingridxlan/lib/petsc-3.11.3)
set(PETSC_ARCH arch-linux2-cxx-opt)

set(METIS_DIR /home/ingridxlan/lib/metis-5.0.3)

set(HDF5_ROOT /home/ingridxlan/lib/hdf5-1.8.16)

set(BESSEL_DIR /home/ingridxlan/lib/complex_bessel)
MESSAGE(STATUS "Bessel fun Dir = " ${BESSEL_DIR})

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})

include_directories(${BESSEL_DIR}/include)
include_directories(${BESSEL_DIR}/include/complex_bessel_bits)

LINK_DIRECTORIES( ${BESSEL_DIR}/lib )
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} complex_bessel)

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

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  /home/ingridxlan/lib/petsc-3.11.3/arch-linux2-cxx-opt/bin/mpicc)
set(CMAKE_CXX_COMPILER /home/ingridxlan/lib/petsc-3.11.3/arch-linux2-cxx-opt/bin/mpicxx)
set(CMAKE_CXX_STANDARD 11)

if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O3 -Wall")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O0 -Wall")
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
