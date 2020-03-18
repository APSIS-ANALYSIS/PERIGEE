# THIS IS THE CMAKE FILE THAT DEFINES THE BASIC
# VAIRABLES FOR BUILDING PROJECTS ON THIS SPECIFIC
# MACHINE WITH PETSC AND VTK LIBRARIES.
# NOTE: THIS FILE DEPENDS ON THE MACHINE. THE DEVELOPER
# IS RESPONSIBLE FOR SETTING CORRECT VALUES FOR THESE
# VARIABLES.
# IN THE CMAKELISTS.TXT FILE, YOU ONLY NEED TO INCLUDE
# THIS FILE TO HAVE THESE VARIABLES DEFINED IN YOUR CMAKE.

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# 1. VTK
set(VTK_DIR /home/jliu/lib/VTK-8.1.1/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})

# 2. PETSc & METIS
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-opt)
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set(PETSC_DIR /home/jliu/lib/petsc-3.9.3-debug)
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

set(PETSC_ARCH .)

set(METIS_DIR /home/jliu/lib/metis-5.0.3)

# ---------------------------------------------------------
find_package(PETSc)

include_directories(${PETSC_INC})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})

if(PETSC_METIS)
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
else(PETSC_METIS)
  find_package(METIS)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LIBRARIES})
endif(PETSC_METIS)

# 3. HDF5
SET(HDF5_ROOT /home/jliu/lib/hdf5-1.8.20)

find_package( HDF5 )

include_directories( ${HDF5_INCLUDE_DIRS} )

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

# ========================================================
# 5. EIGEN variables 
# ========================================================
#INCLUDE_DIRECTORIES(/home/jliu/lib/eigen-3.3.4/include/eigen3)

# ========================================================
# 6. Slepc VARIABLES
# ========================================================
IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-opt)
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(SLEPC_DIR /home/jliu/lib/slepc-3.9.2-debug)
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# ========================================================
# 7. Compiler options 
# ========================================================
SET(CMAKE_C_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicc)
SET(CMAKE_CXX_COMPILER  /home/jliu/lib/mpich-3.2.1/bin/mpicxx)
IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O3 -Wall")
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-DENABLE_TEST -O0 -Wall")
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
SET(CMAKE_CXX_STANDARD 11)

message(STATUS ${EXTRA_LINK_LIBS})

# This an example for buiding code with prof for profiling
#SET(CMAKE_CXX_FLAGS "-pg -Wall")
