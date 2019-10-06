# =========================================================
# 1. VTK VARIABLES
# =========================================================
SET(VTK_DIR /home/groups/amarsden/lib-perigee/VTK-7.1.1-shared)

#----------------------------------------------------------
SET(VTK_VERSION vtk-7.1)
SET(VTK_link_lib vtkCommonCore-7.1 vtkCommonSystem-7.1 vtkCommonDataModel-7.1
  vtkCommonExecutionModel-7.1 vtkCommonMisc-7.1 vtkCommonTransforms-7.1
  vtkCommonMath-7.1 vtkIOCore-7.1 vtkIOLegacy-7.1 vtkIOXML-7.1 vtksys-7.1 
  vtkzlib-7.1 )
#----------------------------------------------------------

# ========================================================
# 2. PETSc VARIABLES
# ========================================================
SET(PETSC_DIR /home/groups/amarsden/lib-perigee/petsc-3.11.3)
SET(PETSC_ARCH arch-linux2-c-opt)

#SET(PETSC_DIR /home/groups/amarsden/lib-perigee/petsc-3.8.4-opt)
#SET(PETSC_ARCH .)

#----------------------------------------------------------
# You do not need to modify the rest of PETSc variables.
SET(PETSC_LIBRARY_DIRS ${PETSC_DIR}/${PETSC_ARCH}/lib )

find_library (PETSC_LIBRARIES NAMES petsc HINTS "${PETSC_DIR}/${PETSC_ARCH}" 
  PATH_SUFFIXES "lib" NO_DEFAULT_PATH)

find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}/${PETSC_ARCH}"
  PATH_SUFFIXES "lib/petsc/conf" "conf" NO_DEFAULT_PATH)

include(${PETSC_CONF_DIR}/PETScBuildInternal.cmake)
#include(${PETSC_CONF_DIR}/PETScConfig.cmake)

SET(PETSC_link_lib ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS})
#----------------------------------------------------------

# ========================================================
# 3. METIS VARIABLES
# ========================================================
SET(METIS_DIR /home/groups/amarsden/lib-perigee/metis-5.0.3)

# ========================================================
# 4. HDF5 VARIABLES
# ========================================================
# Modify the HDF5_DIR
SET(HDF5_DIR /home/groups/amarsden/lib-perigee/hdf5-1.8.16)

# ========================================================
# 5. Compiler options 
# ========================================================
SET(CMAKE_C_COMPILER  /home/groups/amarsden/lib-perigee/petsc-3.11.3/arch-linux2-c-opt/bin/mpicc)
SET(CMAKE_CXX_COMPILER /home/groups/amarsden/lib-perigee/petsc-3.11.3/arch-linux2-c-opt/bin/mpicxx)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O3 -Wall")
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O0 -Wall")
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# EOF
