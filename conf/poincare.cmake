# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# This one is for poincare, my laptop.
# =========================================================
# 1. VTK VARIABLES
# =========================================================
SET(VTK_DIR /Users/juliu/lib/VTK-7.1.1-debug)
SET(VTK_VERSION vtk-7.1)
SET(VTK_link_lib vtkCommonCore-7.1 vtkCommonSystem-7.1 vtkCommonDataModel-7.1
  vtkCommonExecutionModel-7.1 vtkCommonMisc-7.1 vtkCommonTransforms-7.1
  vtkCommonMath-7.1 vtkIOCore-7.1 vtkIOLegacy-7.1 vtkIOXML-7.1 vtksys-7.1
  vtkzlib-7.1 )

# ========================================================
# 2. PETSc VARIABLES
# ========================================================
set(PETSC_DIR /Users/juliu/lib/petsc-3.11.3-debug)

set(PETSC_ARCH .)

find_package(PETSc)

#SET(PETSC_LIBRARY_DIRS ${PETSC_DIR}/${PETSC_ARCH}/lib )
#find_library (PETSC_LIBRARIES NAMES petsc HINTS "${PETSC_DIR}" 
#  PATH_SUFFIXES "${PETSC_ARCH}/lib" "lib" NO_DEFAULT_PATH)
#find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}/${PETSC_ARCH}"
#  PATH_SUFFIXES "lib/petsc/conf" "conf" NO_DEFAULT_PATH)
#include(${PETSC_CONF_DIR}/PETScBuildInternal.cmake)
#include(${PETSC_CONF_DIR}/PETScConfig.cmake)
#SET(PETSC_link_lib ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS})


# ========================================================
# 3. METIS VARIABLES
# ========================================================
#SET(METIS_DIR ${PETSC_DIR})

# ========================================================
# 4. HDF5 VARIABLES
# ========================================================
SET(HDF5_DIR ${PETSC_DIR})

# ========================================================
# 5. Eigen VARIABLES
# ========================================================
#INCLUDE_DIRECTORIES(/Users/juliu/lib/petsc-3.9.4-debug/include/eigen3)

# ========================================================
# 5. Compiler options 
# ========================================================
SET(CMAKE_C_COMPILER  /Users/juliu/lib/mpich-3.2.1/bin/mpicc)
SET(CMAKE_CXX_COMPILER  /Users/juliu/lib/mpich-3.2.1/bin/mpicxx)
SET(CMAKE_CXX_FLAGS "-O0 -Wall")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
