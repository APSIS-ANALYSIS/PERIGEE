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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

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

set(PETSC_DIR /work/01346/liujuy/stampede2/lib/petsc-3.6.4-opt)

set(PETSC_ARCH .)

find_package(PETSc)

# ===================================================================
# 4. HDF5 VARIABLES
# ===================================================================
SET(HDF5_DIR /opt/apps/intel17/impi17_0/phdf5/1.8.16/x86_64 ) 

# ===================================================================
# 5. Compiler options
# ===================================================================
SET(CMAKE_C_COMPILER  /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicc)
SET(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpicxx)
SET(CMAKE_CXX_FLAGS "-xCORE-AVX512 -O3 -xhost -Wall")
SET(CMAKE_BUILD_TYPE RELEASE)

# ===================================================================
# END OF FILE
# ===================================================================
