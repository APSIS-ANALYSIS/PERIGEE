project(perigee)

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

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
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/System/Matrix_3x3.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/Gmsh_FileIO.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P1.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P2.cpp
  ${perigee_source}/Mesh/INodalBC.cpp
  ${perigee_source}/Mesh/IElemBC.cpp
  ${perigee_source}/Mesh/Mesh_Tet4.cpp
  ${perigee_source}/Mesh/Mesh_Tet10.cpp
  ${perigee_source}/Mesh/IPart.cpp
  ${perigee_source}/Mesh/Part_Tet.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet_wall.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet_outflow.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D_inflow.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp_wall.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp_outflow.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  )

SET( perigee_analysis_lib_src
  ${perigee_source}/System/PETSc_Tools.cpp
  ${perigee_source}/System/Matrix_double_6by6_Array.cpp
  ${perigee_source}/System/Matrix_PETSc.cpp
  ${perigee_source}/System/Tensor4_3D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Mesh/QuadPts_debug.cpp
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_FEM_3D.cpp
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Inflow_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC_outflow.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC_wall.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Element/FEANode.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  ${perigee_source}/Element/FEAElement_Tet10_v2.cpp
  ${perigee_source}/Element/FEAElement_Triangle6.cpp
  ${perigee_source}/Element/FEAElement_Tet10.cpp
  ${perigee_source}/Element/FEAElement_Triangle3.cpp
  ${perigee_source}/Model/GenBC_RCR.cpp
  ${perigee_source}/Model/GenBC_Resistance.cpp
  ${perigee_source}/Model/GenBC_Inductance.cpp
  ${perigee_source}/Model/CVFlowRate_Unsteady.cpp
  ${perigee_source}/Model/CVFlowRate_Linear2Steady.cpp
  ${perigee_source}/Model/MaterialModel_StVenant_Kirchhoff.cpp
  ${perigee_source}/Model/MaterialModel_StVenant_Kirchhoff_Simo85.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_RWR2000.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_ST91.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_Incompressible_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_StVenant_Kirchhoff_M94_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_M94_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_Quad_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_GOH06_ST91_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_ST91_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_Guccione_Incompressible_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_GOH06_Incompressible_Mixed.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  )

SET( perigee_postprocess_lib_src
  ${perigee_source}/Mesh/QuadPts_vis_tri6.cpp
  ${perigee_source}/Mesh/QuadPts_vis_tet4.cpp
  ${perigee_source}/Mesh/QuadPts_vis_tet10_v2.cpp
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  )

INCLUDE_DIRECTORIES( ../../include )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess ${EXTRA_LINK_LIBS} )

# 2. Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis ${EXTRA_LINK_LIBS} perigee_preprocess )

# 3. Postprocess lib
ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis
  perigee_preprocess )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( sys_test sys_test.cpp )
ADD_EXECUTABLE( mod_test matmodel_test.cpp )
ADD_EXECUTABLE( testh5 h5_test.cpp )
ADD_EXECUTABLE( elem_test elem_test.cpp )

TARGET_LINK_LIBRARIES( sys_test perigee_preprocess )
TARGET_LINK_LIBRARIES( mod_test perigee_preprocess perigee_analysis )
TARGET_LINK_LIBRARIES( testh5 perigee_preprocess )
TARGET_LINK_LIBRARIES( elem_test perigee_analysis )

# EOF
