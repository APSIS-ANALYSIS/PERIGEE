# This is a sample code that provide a linkage between the
# PERIGEE code and external libraries. Users will have to
# provide correct values for the libraries, according to their
# installation.
# This one is a sample one, assuming the libraries are installed
# following the guide documented in
# https://github.com/ju-liu/PERIGEE-NS/blob/master/install-external-libs.md
# The value of $HOME is /home/jliu

# ========================================================
# Specify the library locations
# ========================================================
# VTK_DIR should be the vtk directory/lib/cmake/vtk-version,
# which contains VTKConfig.cmake file.
# In the guide, vtk directory is $HOME/lib/VTK-7.1.1-shared
set(VTK_DIR /Users/zinanhu/lib/VTK-9.0.1-shared/lib/cmake/vtk-9.0)

# Modify the PETSC_DIR variable to point to the location of PETSc.
set(PETSC_DIR /Users/zinanhu/lib/petsc)

# Modify the PETSC_ARCH variable. You can find it in your configuration
# output. If you forget it, go to your PETSc home director and open
# configure.log. Go the end of the file, and you shall find the value 
# of PETSC_ARCH
set(PETSC_ARCH arch-darwin-cxx-debug)

# Modify the METIS_DIR.
# Note: If your PETSc has METIS installed, the conf
# file will directly load that METIS; otherwise this METIS will
# be used for PERIGEE. This means, if you are sure that you have
# METIS in PETSc, you do not have to specify the METIS_DIR variable.
set(METIS_DIR /Users/zinanhu/lib/metis-5.1.0)

# Modify the HDF5_ROOT, pointing to your hdf5 library location
set(HDF5_ROOT /Users/zinanhu/lib/hdf5-1.8.16)

# ========================================================
# Setup the libraries
# You do NOT have to modify anything in this part
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

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

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
# Specify the MPI compilers. There should be compilers in
# $PETSC_DIR/$PETSC_ARCH/bin, or the mpich you specified for 
# PETSc install.
set(CMAKE_C_COMPILER  /Users/zinanhu/lib/petsc/arch-darwin-cxx-debug/bin/mpicc)
set(CMAKE_CXX_COMPILER /Users/zinanhu/lib/petsc/arch-darwin-cxx-debug/bin/mpicxx)
set(CMAKE_CXX_STANDARD 11)
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-O3 -Wall")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
