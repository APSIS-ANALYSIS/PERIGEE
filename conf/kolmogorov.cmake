# Configuration setup for machine kolmogorov mac laptop

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR /Users/juliu/lib/VTK-8.2.0/lib/cmake/vtk-8.2)

set(PETSC_DIR /Users/juliu/lib/petsc-3.16.6-debug)

set(PETSC_ARCH .)

set(METIS_DIR /Users/juliu/lib/metis-5.1.0)

set(HDF5_ROOT /Users/juliu/lib/hdf5-1.12.0)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy vtkIOXML REQUIRED)
find_package(HDF5 REQUIRED)
find_package(PETSc REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PETSC_INC})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} -L/Users/juliu/lib/petsc-3.16.6-debug/lib -Wl,-rpath,/Users/juliu/lib/petsc-3.16.6-debug/lib -L/Users/juliu/lib/petsc-3.16.6-debug/lib -Wl,-rpath,/usr/local/lib -L/usr/local/lib -Wl,-rpath,/Users/juliu/lib/mpich-3.4rc1/lib -L/Users/juliu/lib/mpich-3.4rc1/lib -Wl,-rpath,/usr/local/Cellar/gcc/12.2.0/lib/gcc/current/gcc/x86_64-apple-darwin19/12 -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/current/gcc/x86_64-apple-darwin19/12 -Wl,-rpath,/usr/local/Cellar/gcc/12.2.0/lib/gcc/current/gcc -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/current/gcc -Wl,-rpath,/usr/local/Cellar/gcc/12.2.0/lib/gcc/current -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/current -lpetsc -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lchaco -lparmetis -lmetis -lyaml -lX11 -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lgcc_s.1.1 -lquadmath -lc++ -ldl)

if(PETSC_METIS)
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
  message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
else(PETSC_METIS)
  find_package(METIS)
  include_directories(${METIS_INCLUDE_DIRS})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LIBRARIES})
endif(PETSC_METIS)

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
set(CMAKE_C_COMPILER  /Users/juliu/lib/mpich-3.4rc1/bin/mpicc)
set(CMAKE_CXX_COMPILER  /Users/juliu/lib/mpich-3.4rc1/bin/mpicxx)
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
