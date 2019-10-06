# This is a sample code that provide a linkage between the
# PERIGEE code and external libraries. Users will have to
# provide correct values for the libraries, according to their
# installation.
# This one is a sample one, assuming the libraries are installed
# following the guide documented in
# https://github.com/ju-liu/PERIGEE-NS/blob/master/install-external-libs.md
# The value of $HOME will be /home/jliu

# =========================================================
# 1. VTK VARIABLES
# =========================================================
# VTK_DIR should be the prefix value you run cmake when installing VTK.
# In the guide, it is $HOME/lib/VTK-7.1.1-shared
SET(VTK_DIR /home/ingridxlan/lib/VTK-7.1.1-shared)

# You do not need to modify VTK_VERSION and VTK_link_lib, if you are 
# using VTK-7.x.x
SET(VTK_VERSION vtk-7.1)
SET(VTK_link_lib vtkCommonCore-7.1 vtkCommonSystem-7.1 vtkCommonDataModel-7.1
  vtkCommonExecutionModel-7.1 vtkCommonMisc-7.1 vtkCommonTransforms-7.1
  vtkCommonMath-7.1 vtkIOCore-7.1 vtkIOLegacy-7.1 vtkIOXML-7.1 vtksys-7.1 
  vtkzlib-7.1 )

# ========================================================
# 2. PETSc VARIABLES
# ========================================================
# Modify the PETSC_DIR variable to point to the location of PETSc.
SET(PETSC_DIR /home/ingridxlan/lib/petsc-3.11.3)

# Modify the PETSC_ARCH variable. You can find it in your configuration
# output. If you forget it, go to your PETSc home director and open
# configure.log. Go the end of the file, and you shall find the value 
# of PETSC_ARCH
SET(PETSC_ARCH arch-linux2-c-debug)

# You do not need to modify the rest of PETSc variables.
SET(PETSC_LIBRARY_DIRS ${PETSC_DIR}/${PETSC_ARCH}/lib )

find_library (PETSC_LIBRARIES NAMES petsc HINTS "${PETSC_DIR}/${PETSC_ARCH}" 
  PATH_SUFFIXES "lib" NO_DEFAULT_PATH)

find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}/${PETSC_ARCH}"
  PATH_SUFFIXES "lib/petsc/conf" "conf" NO_DEFAULT_PATH)

include(${PETSC_CONF_DIR}/PETScBuildInternal.cmake)
#include(${PETSC_CONF_DIR}/PETScConfig.cmake)

SET(PETSC_link_lib ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS})

# ========================================================
# 3. METIS VARIABLES
# ========================================================
# Modify the METIS_DIR.
SET(METIS_DIR /home/ingridxlan/lib/metis-5.0.3)

# ========================================================
# 4. HDF5 VARIABLES
# ========================================================
# Modify the HDF5_DIR
SET(HDF5_DIR /home/ingridxlan/lib/hdf5-1.8.16)

# ========================================================
# 5. Compiler options 
# ========================================================
# Specify the MPI compilers. There should be compilers in
# $PETSC_DIR/$PETSC_ARCH/bin
SET(CMAKE_C_COMPILER  /home/ingridxlan/lib/petsc-3.11.3/arch-linux2-c-debug/bin/mpicc)
SET(CMAKE_CXX_COMPILER /home/ingridxlan/lib/petsc-3.11.3/arch-linux2-c-debug/bin/mpicxx)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O3 -Wall")
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O0 -Wall")
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# EOF
