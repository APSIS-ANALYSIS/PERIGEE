# ===================================================================
# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.
# -------------------------------------------------------------------
# This one is for stampede2 build.
# ===================================================================

# ===================================================================
# 1. VTK VARIABLES
# ===================================================================
SET(VTK_DIR /work/01346/liujuy/stampede2/lib/VTK-7.1.1)
SET(VTK_VERSION vtk-7.1)
SET(VTK_link_lib vtkCommonCore-7.1 vtkCommonSystem-7.1 vtkIOCore-7.1 vtksys-7.1 vtkCommonDataModel-7.1 vtkIOXML-7.1 vtkIOLegacy-7.1 vtkCommonExecutionModel-7.1 vtkCommonMisc-7.1 vtkCommonTransforms-7.1 vtkCommonMath-7.1 vtkIOCore-7.1 vtkzlib-7.1 )

MESSAGE(STATUS ${VTK_DIR})

# ===================================================================
# 2. PETSc VARIABLES
# ===================================================================
#SET(PETSC_DIR /home1/apps/intel18/impi18_0/petsc/3.10)
#SET(PETSC_ARCH skylake)

SET(PETSC_DIR /work/01346/liujuy/stampede2/lib/petsc-3.6.4-opt)
SET(PETSC_ARCH .)

SET(PETSC_LIBRARY_DIRS ${PETSC_DIR}/${PETSC_ARCH}/lib )

find_library (PETSC_LIBRARIES NAMES petsc HINTS "${PETSC_DIR}/${PETSC_ARCH}" 
  PATH_SUFFIXES "lib" NO_DEFAULT_PATH)

find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}/${PETSC_ARCH}"
  PATH_SUFFIXES "lib/petsc/conf" "conf" NO_DEFAULT_PATH)

include(${PETSC_CONF_DIR}/PETScConfig.cmake)
#include(${PETSC_CONF_DIR}/PETScBuildInternal.cmake)

# Load libpetsc and the external libs
SET(PETSC_link_lib ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS})

# ===================================================================
# 3. METIS VARIABLES
# ===================================================================
SET(METIS_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib)

# ===================================================================
# 4. HDF5 VARIABLES
# ===================================================================
SET(HDF5_DIR /opt/apps/intel17/impi17_0/phdf5/1.8.16/x86_64 ) 

# ===================================================================
# 5. Eigen VARIABLES
# ===================================================================
#INCLUDE_DIRECTORIES(/work/01346/liujuy/stampede2/lib/eigen-3.3.4/include/eigen3)

# ===================================================================
# 6. SLEPc VARIABLES
# ===================================================================
#SET(SLEPC_DIR /work/01346/liujuy/stampede2/lib/slepc-3.8.3-opt)

# ===================================================================
# 7. Compiler options
# ===================================================================
SET(CMAKE_C_COMPILER  /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicc)
SET(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicxx)
SET(CMAKE_CXX_FLAGS "-xCORE-AVX512 -O3 -xhost -Wall")
SET(CMAKE_BUILD_TYPE RELEASE)

# ===================================================================
# END OF FILE
# ===================================================================
