PROJECT(perigee)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
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

# print makefile details or not
SET(CMAKE_VERBOSE_MAKEFILE FALSE)

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)
INCLUDE(../../conf/setup_external_libs.cmake)

SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vec_Tools.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/Gmsh_FileIO.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P1.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P2.cpp
  ${perigee_source}/Mesh/INodalBC.cpp
  ${perigee_source}/Mesh/IElemBC.cpp
  ${perigee_source}/Mesh/Mesh_Tet4.cpp
  ${perigee_source}/Mesh/Mesh_Tet10.cpp
  ${perigee_source}/Mesh/IPart.cpp
  ${perigee_source}/Mesh/Part_Tet4.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet4.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet4_outflow.cpp
  ${perigee_source}/Mesh/NodalBC_3D_vtp.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D_inflow.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp_outflow.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  )


INCLUDE_DIRECTORIES( ../../include )
INCLUDE_DIRECTORIES( ./include )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess ${EXTRA_LINK_LIBS} )

# 2. Analysis libs
#ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
#TARGET_LINK_LIBRARIES( perigee_analysis ${EXTRA_LINK_LIBS} )

# 3. Postprocess lib
#ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
#TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)




ADD_EXECUTABLE( preprocess3d preprocess_sv_tets.cpp)

TARGET_LINK_LIBRARIES( preprocess3d ${EXTRA_LINK_LIBS} perigee_preprocess )

ADD_EXECUTABLE( testtet test_tet_t.cpp)

TARGET_LINK_LIBRARIES( testtet ${EXTRA_LINK_LIBS} perigee_preprocess )

ADD_EXECUTABLE( testtetsurface test_tet_t_vtp.cpp)

TARGET_LINK_LIBRARIES( testtetsurface ${EXTRA_LINK_LIBS} perigee_preprocess )

ADD_EXECUTABLE( vtu2vtp vtu2vtp.cpp)
TARGET_LINK_LIBRARIES( vtu2vtp ${EXTRA_LINK_LIBS} perigee_preprocess )

