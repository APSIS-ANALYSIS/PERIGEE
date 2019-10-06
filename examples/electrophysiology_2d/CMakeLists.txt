PROJECT(perigee)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

SET(perigee_VERSION_MAJOR 0)
SET(perigee_VERSION_MINOR 0)

# print on screen the system type
MESSAGE(STATUS "THE SYSTEM is " ${CMAKE_SYSTEM_NAME})
MESSAGE(STATUS "              " ${CMAKE_SYSTEM})

SET(CMAKE_BUILD_TYPE Debug CACHE STRING "default to debug" FORCE)
OPTION(BUILD_SHARED_LIBS "build shard libraries" OFF)

# Libraries location
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# binary location
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# print makefile details or not
SET(CMAKE_VERBOSE_MAKEFILE FALSE)

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)

MESSAGE(STATUS "=================== External library setup ===================")
# ===============================================
# SET VTK
# ===============================================
IF(DEFINED VTK_DIR)
  MESSAGE(STATUS "VTK_DIR = ${VTK_DIR}")
ELSE(DEFINED VTK_DIR)
  MESSAGE(FATAL_ERROR "VTK_DIR is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_DIR)

If(DEFINED VTK_VERSION)
  MESSAGE(STATUS "VTK_VERSION = ${VTK_VERSION}")
ELSE(DEFINED VTK_VERSION)
  MESSAGE(FATAL_ERROR "VTK_VERSION is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_VERSION)

IF(DEFINED VTK_link_lib)
  MESSAGE(STATUS "VTK_link_lib=${VTK_link_lib}")
ELSE(DEFINED VTK_link_lib)
  MESSAGE(FATAL_ERROR "VTK_link_lib is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_link_lib)

SET(VTK_INCLUDE_DIRS ${VTK_INCLUDE_DIRS}
  ${VTK_DIR}/include/${VTK_VERSION})

IF(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib/${VTK_VERSION})
ELSE(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib)
ENDIF(${VTK_VERSION} MATCHES "vtk-5.10")

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_link_lib} )
# ===============================================
# FINISH VTK SETUP
# ===============================================

# ===============================================
# SET PETSc library 
# ===============================================
IF(DEFINED PETSC_DIR)
  MESSAGE(STATUS "PETSc directory = ${PETSC_DIR}")
ELSE(DEFINED PETSC_DIR)
  MESSAGE(FATAL_ERROR "PETSC_DIR is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_DIR)

IF(DEFINED PETSC_ARCH)
  MESSAGE(STATUS "PETSc arch = ${PETSC_ARCH}")
ELSE(DEFINED PETSC_ARCH)
  MESSAGE(FATAL_ERROR "PETSC_ARCH is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_ARCH)

IF(DEFINED PETSC_LIBRARY_DIRS)
  MESSAGE(STATUS "PETSC_LIBRARY_DIRS = ${PETSC_LIBRARY_DIRS}")
ELSE(DEFINED PETSC_LIBRARY_DIRS)
  MESSAGE(FATAL_ERROR "PETSC_LIBRARY_DIRS is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_LIBRARY_DIRS)

IF(DEFINED PETSC_link_lib)
  MESSAGE(STATUS "PETSC_link_lib = ${PETSC_link_lib}")
ELSE(DEFINED PETSC_link_lib)
  MESSAGE(FATAL_ERROR "PETSC_link_lib is NOT defined! Check
  ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_link_lib)

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_link_lib})

SET(PETSC_INCLUDE_DIRS ${PETSC_DIR}/include
  ${PETSC_DIR}/${PETSC_ARCH}/include)
# ===============================================
# Finish PETSc Library Setup
# ===============================================

# ===============================================
# SET METIS library
# ===============================================
IF(DEFINED METIS_DIR)
  MESSAGE(STATUS "METIS DIR = ${METIS_DIR}")
ELSE(DEFINED METIS_DIR)
  MESSAGE(FATAL_ERROR "METIS_DIR is NOT defined! check ~/PETSc_VTK.cmake")
ENDIF(DEFINED METIS_DIR)

SET(METIS_LIBRARY_DIRS ${METIS_DIR}/lib)
SET(METIS_INCLUDE_DIRS ${METIS_DIR}/include)
SET(METIS_LINK_LIBS metis)

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LINK_LIBS})
# ===============================================
# FINISH METIS setup
# ===============================================

# ===============================================
# SET HDF5 library
# ===============================================
IF(DEFINED HDF5_DIR)
  MESSAGE(STATUS "HDF5 DIR = ${HDF5_DIR}")
ELSE(DEFINED HDF5_DIR)
  MESSAGE(FATAL_ERROR "HDF5_DIR is NOT defined! check ~/PETSc_VTK.cmake")
ENDIF(DEFINED HDF5_DIR)

SET(HDF5_LIBRARY_DIRS ${HDF5_DIR}/lib)
SET(HDF5_INCLUDE_DIRS ${HDF5_DIR}/include)
SET(HDF5_LINK_LIBS hdf5)

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LINK_LIBS})
# ===============================================
# FINISH HDF5 setup
# ===============================================

# ===============================================
# Print Compiler/Linker and flags
# ===============================================
MESSAGE(STATUS "=================== Compiler setup ==========================")
MESSAGE(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER})
MESSAGE(STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS})
MESSAGE(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCC: " ${CMAKE_COMPILER_IS_GNUCC})
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCXX: " ${CMAKE_COMPILER_IS_GNUCXX})
# ===============================================
# FINISH COMPILER SETUP
# ===============================================

# ===============================================
# Set my own source code
# ===============================================
SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/NURBS_Tool/NURBS_Bezier.cpp
  ${perigee_source}/NURBS_Tool/NURBS_Tools.cpp
  ${perigee_source}/NURBS_Tool/NURBS_FileIO.cpp
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vec_Tools.cpp
  ${perigee_source}/System/kRefinement.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/IPart.cpp
  ${perigee_source}/Mesh/Mesh_NURBS_1Patch_2D.cpp
  ${perigee_source}/Mesh/IEN_NURBS_1Patch_2D.cpp
  ${perigee_source}/Mesh/BoundaryCond2D.cpp
  ${perigee_source}/Mesh/BC_Partition2D.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Mesh/Part_NURBS_1Patch_2D_METIS.cpp
  )

# 2. Anlysis source cpp
SET( perigee_analysis_lib_src
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vec_Tools.cpp
  ${perigee_source}/System/Matrix_double_6by6_Array.cpp
  ${perigee_source}/System/Matrix_PETSc.cpp
  ${perigee_source}/System/HDF5_PartReader.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss.cpp
  ${perigee_source}/Mesh/QuadPts_vis.cpp
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_1Patch_NURBS_2D.cpp
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_BC_3D.cpp
  ${perigee_source}/Analysis_Tool/ALocal_BC_2D.cpp
  ${perigee_source}/Analysis_Tool/AInt_Weight.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/ALocal_meshSize_2D_NURBS.cpp
  ${perigee_source}/Element/BernsteinBasis_Array.cpp
  ${perigee_source}/Element/AExtractor_2D_NURBS_xy.cpp
  ${perigee_source}/Element/FEANode.cpp
  ${perigee_source}/Element/FEAElement_NURBS_3D_der0_v3.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der1.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der2.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der0.cpp
  ${perigee_source}/Model/PLocAssem_NLHeat_2D_GenAlpha.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/PDNSolution_heatEqn.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PGAssem.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_source}/Solver/PNonlinear_Solver.cpp
  ${perigee_source}/Solver/PTime_Solver.cpp
  )

# 3. Postprocess cpp
SET( perigee_postprocess_lib_src
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/VTK_Writer.cpp
  )

#  ${perigee_SOURCE_DIR}/src/Post_error.cpp
#  ${perigee_SOURCE_DIR}/src/VisDataPrep_2DNLHeat.cpp

# ===============================================
# Finish setting my own source code
# ===============================================
# INCLUDE LIBRARY HEADER FILES
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${PETSC_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${METIS_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
# LINK TO LIBRARY LIBS
LINK_DIRECTORIES(${VTK_LINK_DIRS})
LINK_DIRECTORIES(${PETSC_LIBRARY_DIRS})
LINK_DIRECTORIES(${METIS_LIBRARY_DIRS})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})

# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )

MESSAGE(STATUS "=================== Misc setup ==========================")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess ${EXTRA_LINK_LIBS} )

# 2. Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis ${EXTRA_LINK_LIBS} )

# 3. Postprocess lib
ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis )
# -------------------------------------------------------------------

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
# 1. Make executable
# 1.1 Preprocess
ADD_EXECUTABLE( preprocess2d preprocess2d_main.cpp)

# 1.2 Analysis and tests
#ADD_EXECUTABLE( analysis_heat heat_2d_genAlpha.cpp)

# 1.3 Post-process
#ADD_EXECUTABLE( prepostproc Pre_postprocess2d.cpp)
#ADD_EXECUTABLE( post_compare_manu post_compare_manu_2d.cpp)
#ADD_EXECUTABLE( vis_2d_nlheat vis_2DNLHeat.cpp)

# -------------------------------------------------------------------
# 2. Link libraries
# 2.1 Preprocess
TARGET_LINK_LIBRARIES( preprocess2d ${EXTRA_LINK_LIBS} 
  perigee_preprocess )

# 2.2 Analysis and tests
#TARGET_LINK_LIBRARIES( analysis_heat ${EXTRA_LINK_LIBS} perigee_analysis )

# 2.3 Postprocess
#TARGET_LINK_LIBRARIES( prepostproc ${EXTRA_LINK_LIBS} perigee_preprocess
#  perigee_analysis )

#TARGET_LINK_LIBRARIES( post_compare_manu ${EXTRA_LINK_LIBS} 
#  perigee_preprocess perigee_analysis perigee_postprocess )

#TARGET_LINK_LIBRARIES( vis_2d_nlheat ${EXTRA_LINK_LIBS}
#  perigee_preprocess perigee_analysis perigee_postprocess )

#EOF
