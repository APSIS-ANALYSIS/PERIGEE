# ===================================================================
# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.
# ===================================================================

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# ===================================================================
# 1. VTK VARIABLES
# ===================================================================
SET(VTK_DIR /work/01346/liujuy/lonestar/lib/VTK-6.3.0)
SET(VTK_VERSION vtk-6.3)
SET(VTK_link_lib vtkCommonCore-6.3 vtkCommonSystem-6.3 vtkIOCore-6.3 vtksys-6.3 vtkCommonDataModel-6.3 vtkIOXML-6.3 vtkIOLegacy-6.3 vtkCommonExecutionModel-6.3 vtkCommonMisc-6.3 vtkCommonTransforms-6.3 vtkCommonMath-6.3 vtkIOCore-6.3 vtkzlib-6.3 )

# ===================================================================
# 2. PETSc VARIABLES
# ===================================================================
#SET(PETSC_DIR /work/01346/liujuy/lonestar/lib/petsc-3.8.4-opt-largemem)
#SET(PETSC_ARCH .)

set(PETSC_DIR /opt/apps/intel18/cray_mpich_7_7/petsc/3.11)
set(PETSC_ARCH haswell)

find_package(PETSc)

# ===================================================================
# 4. HDF5 VARIABLES
# ===================================================================
#SET(HDF5_DIR /opt/apps/intel18/hdf5/1.8.16/x86_64)
SET(HDF5_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib )

# ===================================================================
# 5. Eigen VARIABLES
# ===================================================================
#INCLUDE_DIRECTORIES(/work/01346/liujuy/lonestar/lib/eigen-3.3.4/include/eigen3)

# ===================================================================
# 6. SLEPc VARIABLES
# ===================================================================
#SET(SLEPC_DIR /work/01346/liujuy/lonestar/lib/slepc-3.10.1-opt )

#SET(BOOST_DIR /opt/apps/intel18/boost/1.64)

SET(BESSEL_DIR /work/01346/liujuy/lonestar/lib/complex-bessel)

# ===================================================================
# 6. Compiler options
# ===================================================================
SET(CMAKE_C_COMPILER  /opt/apps/intel18/cray_mpich/7.7.3/bin/mpicc)
SET(CMAKE_CXX_COMPILER  /opt/apps/intel18/cray_mpich/7.7.3/bin/mpicxx)
#SET(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/mpicc)
#SET(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/mpicxx)
SET(CMAKE_CXX_FLAGS "-O3 -xhost -Wall")
SET(CMAKE_BUILD_TYPE RELEASE)
SET(CMAKE_CXX_STANDARD 11)

# ===================================================================
# END OF FILE
# ===================================================================
