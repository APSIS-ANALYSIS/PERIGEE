# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.

# This one is for Stanford bacon Linux machine, generated 
# on Mar. 6 2017. This one is for updateing the PETSc to
# version 3.7.5, and try to get the petsc libs in a automatic
# way by reading the PETScBuildInternal.cmake file.
# The PETScBuildInternal.cmake file will give the libraries
# and their paths for linking.

# =========================================================
# 1. VTK VARIABLES
# =========================================================
SET(VTK_DIR /home/jliu/lib/VTK-8.2.0)
SET(VTK_VERSION vtk-8.2)
SET(VTK_link_lib vtkCommonCore-8.2 vtkCommonSystem-8.2 vtkCommonDataModel-8.2
  vtkCommonExecutionModel-8.2 vtkCommonMisc-8.2 vtkCommonTransforms-8.2
  vtkCommonMath-8.2 vtkIOCore-8.2 vtkIOLegacy-8.2 vtkIOXML-8.2 vtksys-8.2 
  vtkzlib-8.2 vtkFiltersGeometry-8.2 )

# ========================================================
# 2. PETSc VARIABLES
# ========================================================

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(PETSC_DIR /home/jliu/lib/petsc-3.8.4-opt)
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(PETSC_DIR /home/jliu/lib/petsc-3.8.4-debug)
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

SET(PETSC_ARCH .)

# I do not need this variable if I use the PETScBuildInternal.cmake
# file. But my previous loading cmake file requires this variable.
SET(PETSC_LIBRARY_DIRS ${PETSC_DIR}/${PETSC_ARCH}/lib )

# Locate libpetsc
find_library (PETSC_LIBRARIES NAMES petsc HINTS "${PETSC_DIR}/${PETSC_ARCH}" 
  PATH_SUFFIXES "lib" NO_DEFAULT_PATH)

# Locate configuration folder
find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}/${PETSC_ARCH}"
  PATH_SUFFIXES "lib/petsc/conf" "conf" NO_DEFAULT_PATH)

# Load the configuration file
include(${PETSC_CONF_DIR}/PETScBuildInternal.cmake)
#include(${PETSC_CONF_DIR}/PETScConfig.cmake)

# Load libpetsc and the external libs
SET(PETSC_link_lib ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS})

# ========================================================
# 3. METIS VARIABLES
# ========================================================
SET(METIS_DIR /home/jliu/lib/metis-5.0.3-32bit)

# ========================================================
# 4. HDF5 VARIABLES
# ========================================================
SET(HDF5_DIR /home/jliu/lib/hdf5-1.8.17)

# ========================================================
# 5. Eigen VARIABLES
# ========================================================
INCLUDE_DIRECTORIES(/home/jliu/lib/Eigen-3.3.7/include/eigen3)

# ========================================================
# 6. Slepc VARIABLES
# ========================================================
#IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
#  SET(SLEPC_DIR /home/jliu/lib/slepc-3.10.2-opt)
#ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
#  SET(SLEPC_DIR /home/jliu/lib/slepc-3.10.2-debug)
#ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# ========================================================
# 7. Compiler options 
# ========================================================
SET(CMAKE_C_COMPILER  /home/jliu/lib/mpich-3.3.2/bin/mpicc)
SET(CMAKE_CXX_COMPILER /home/jliu/lib/mpich-3.3.2/bin/mpicxx)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O3 -Wall")
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# EOF
