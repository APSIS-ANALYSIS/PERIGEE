CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

project(perigee)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0074 NEW)
ENDIF(COMMAND cmake_policy)

SET(perigee_VERSION_MAJOR 0)
SET(perigee_VERSION_MINOR 0)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "set to release" FORCE)
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "default to debug" FORCE)
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

OPTION(BUILD_SHARED_LIBS "build shard libraries" OFF)

# Libraries location
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# binary location
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)

SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
# 1. Preprocessor source cpp
SET( perigee_analysis_lib_src
  ${perigee_source}/System/PETSc_Tools.cpp
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/Tensor2_3D.cpp
  ${perigee_source}/System/Tensor4_3D.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_1D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Hex.cpp
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_FEM_3D.cpp
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/FEANode.cpp
  ${perigee_source}/Element/FE_Tools.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  ${perigee_source}/Element/FEAElement_Tet10_v2.cpp
  ${perigee_source}/Element/FEAElement_Hex8.cpp
  ${perigee_source}/Element/FEAElement_Hex27.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_Stress_Recovery.cpp
  ${perigee_SOURCE_DIR}/src/PGAssem_Stress_Recovery.cpp
  )

SET( perigee_postprocess_lib_src
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_Stress_Recovery.cpp
  ${perigee_SOURCE_DIR}/src/VTK_Writer_Stress_Recovery.cpp
  )

INCLUDE_DIRECTORIES( ../../include )
INCLUDE_DIRECTORIES( ./include )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis ${EXTRA_LINK_LIBS} )

ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( smooth driver.cpp)
ADD_EXECUTABLE( vis_smooth vis_smooth.cpp )

TARGET_LINK_LIBRARIES( smooth perigee_analysis )
TARGET_LINK_LIBRARIES( vis_smooth perigee_postprocess )

# EOF
