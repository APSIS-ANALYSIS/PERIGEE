PROJECT(perigee)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

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

# Libraries location
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# binary location
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)

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
  ${perigee_source}/Mesh/Mesh_NURBS_1Patch_2D.cpp
  ${perigee_source}/Mesh/Mesh_NURBS_1Patch_3D.cpp
  ${perigee_source}/Mesh/Mesh_NURBS_1ofnPatch_3D.cpp
  ${perigee_source}/Mesh/IEN_NURBS_1Patch_3D.cpp
  ${perigee_source}/Mesh/IEN_NURBS_1Patch_2D.cpp
  ${perigee_source}/Mesh/BoundaryCond2D.cpp
  ${perigee_source}/Mesh/BoundaryCond.cpp
  ${perigee_source}/Mesh/BC_Partition2D.cpp
  ${perigee_source}/Mesh/BC_Partition.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Mesh/Part_NURBS_1Patch_2D_METIS.cpp
  ${perigee_source}/Mesh/Part_NURBS_1Patch_3D_METIS.cpp
  ${perigee_source}/Mesh/IPart.cpp
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
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_1Patch_NURBS_3D.cpp
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_BC_3D.cpp
  ${perigee_source}/Analysis_Tool/ALocal_BC_2D.cpp
  ${perigee_source}/Analysis_Tool/AInt_Weight.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/ALocal_meshSize_2D_NURBS.cpp
  ${perigee_source}/Analysis_Tool/ALocal_meshSize_3D_NURBS.cpp
  ${perigee_source}/Element/BernsteinBasis_Array.cpp
  ${perigee_source}/Element/AExtractor_2D_NURBS_xy.cpp
  ${perigee_source}/Element/AExtractor_3D_NURBS_xyz.cpp
  ${perigee_source}/Element/FEANode.cpp
  ${perigee_source}/Element/FEAElement_NURBS_3D_der0_v3.cpp
  ${perigee_source}/Element/FEAElement_NURBS_3D_der1_v3.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der1.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der2.cpp
  ${perigee_source}/Element/FEAElement_NURBS_2D_der0.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PGAssem.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_source}/Solver/PNonlinear_Solver.cpp
  ${perigee_source}/Solver/PTime_Solver.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_NLHeat_2D_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_NLHeat_3D_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/src/PDNSolution_heatEqn.cpp
  )

# 3. Postprocess cpp
SET( perigee_postprocess_lib_src
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/VTK_Writer.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_2DNLHeat.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_3DNLHeat.cpp
  ${perigee_SOURCE_DIR}/src/Post_error.cpp
  )
# ===============================================
# Finish setting my own source code
# ===============================================

# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )
INCLUDE_DIRECTORIES( ./include )

MESSAGE(STATUS "=================== Misc setup ==========================")
MESSAGE(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
MESSAGE(STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
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
ADD_EXECUTABLE( preprocess3d preprocess3d_main.cpp)

# 1.2 Analysis and tests
ADD_EXECUTABLE( analysis_2d_heat heat_2d_genAlpha.cpp)
ADD_EXECUTABLE( analysis_3d_heat heat_3d_genAlpha.cpp)

# 1.3 Post-process
ADD_EXECUTABLE( prepostproc2d Pre_postprocess2d.cpp)
ADD_EXECUTABLE( post_compare_manu_2d post_compare_manu_2d.cpp)
ADD_EXECUTABLE( vis_2d_nlheat vis_2DNLHeat.cpp)

ADD_EXECUTABLE( prepostproc3d Pre_postprocess3d.cpp)
ADD_EXECUTABLE( post_compare_manu_3d post_compare_manu_3d.cpp)
ADD_EXECUTABLE( vis_3d_nlheat vis_3DNLHeat.cpp)

# -------------------------------------------------------------------
# 2. Link libraries
# 2.1 Preprocess
TARGET_LINK_LIBRARIES( preprocess2d ${EXTRA_LINK_LIBS} perigee_preprocess )
TARGET_LINK_LIBRARIES( preprocess3d ${EXTRA_LINK_LIBS} perigee_preprocess )

# 2.2 Analysis and tests
TARGET_LINK_LIBRARIES( analysis_2d_heat ${EXTRA_LINK_LIBS} perigee_analysis )
TARGET_LINK_LIBRARIES( analysis_3d_heat ${EXTRA_LINK_LIBS} perigee_analysis )

# 2.3 Postprocess
TARGET_LINK_LIBRARIES( prepostproc2d ${EXTRA_LINK_LIBS} perigee_preprocess
  perigee_analysis )
TARGET_LINK_LIBRARIES( prepostproc3d ${EXTRA_LINK_LIBS} perigee_preprocess
  perigee_analysis )

TARGET_LINK_LIBRARIES( post_compare_manu_2d ${EXTRA_LINK_LIBS} 
  perigee_preprocess perigee_analysis perigee_postprocess )

TARGET_LINK_LIBRARIES( post_compare_manu_3d ${EXTRA_LINK_LIBS} 
  perigee_preprocess perigee_analysis perigee_postprocess )

TARGET_LINK_LIBRARIES( vis_2d_nlheat ${EXTRA_LINK_LIBS}
  perigee_preprocess perigee_analysis perigee_postprocess )

TARGET_LINK_LIBRARIES( vis_3d_nlheat ${EXTRA_LINK_LIBS}
  perigee_preprocess perigee_analysis perigee_postprocess )

#EOF
