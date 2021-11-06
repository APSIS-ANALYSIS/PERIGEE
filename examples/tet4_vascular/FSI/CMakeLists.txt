PROJECT(perigee)

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

INCLUDE(../../../conf/basic_variable_log.cmake)
INCLUDE(../../../conf/system_lib_loading.cmake)

# ===============================================
# Set my own source code
# ===============================================
SET( perigee_source ${perigee_SOURCE_DIR}/../../../src )
# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/Matrix_3x3.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/Gmsh_FileIO.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P1.cpp
  ${perigee_source}/Mesh/INodalBC.cpp
  ${perigee_source}/Mesh/IElemBC.cpp
  ${perigee_source}/Mesh/Mesh_Tet4.cpp
  ${perigee_source}/Mesh/IPart.cpp
  ${perigee_source}/Mesh/Part_Tet.cpp
  ${perigee_source}/Mesh/Part_Tet_FSI.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet_outflow.cpp
  ${perigee_source}/Mesh/NodalBC_3D_vtp.cpp
  ${perigee_source}/Mesh/NodalBC_3D_vtu.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/NBC_Partition.cpp
  ${perigee_source}/Mesh/NBC_Partition_inflow.cpp
  ${perigee_source}/Mesh/EBC_Partition.cpp
  ${perigee_source}/Mesh/EBC_Partition_outflow.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  )

# 2. Anlysis source cpp
SET( perigee_analysis_lib_src
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/PETSc_Tools.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/System/Matrix_PETSc.cpp
  ${perigee_source}/System/Matrix_3x3.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/System/Tensor4_3D.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_FEM_3D.cpp
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem_wTag.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Inflow_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC_outflow.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/APart_Node_FSI.cpp
  ${perigee_source}/Element/FEANode.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_M94_Mixed.cpp
  ${perigee_source}/Model/MaterialModel_NeoHookean_Incompressible_Mixed.cpp
  ${perigee_source}/Model/Prestress_solid.cpp
  ${perigee_source}/Model/GenBC_RCR.cpp
  ${perigee_source}/Model/GenBC_Resistance.cpp
  ${perigee_source}/Model/GenBC_Inductance.cpp
  ${perigee_source}/Model/GenBC_Coronary.cpp
  ${perigee_source}/Model/GenBC_Pressure.cpp
  ${perigee_source}/Model/GenBC_Tools.cpp
  ${perigee_source}/Model/CVFlowRate_Steady.cpp
  ${perigee_source}/Model/CVFlowRate_Unsteady.cpp
  ${perigee_source}/Model/CVFlowRate_Linear2Steady.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/Seg_Sol_Tools.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/../src/PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/../src/PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/../src/PLocAssem_Tet4_VMS_Seg_Incompressible.cpp
  ${perigee_SOURCE_DIR}/../src/PLocAssem_Tet4_FSI_Mesh_Elastostatic.cpp
  ${perigee_SOURCE_DIR}/../src/PGAssem_FSI_FEM.cpp
  ${perigee_SOURCE_DIR}/../src/PGAssem_Seg_FEM.cpp
  ${perigee_SOURCE_DIR}/../src/PDNSolution_Mixed_UPV_3D.cpp
  ${perigee_SOURCE_DIR}/../src/PDNSolution_P_V_Mixed_3D.cpp
  ${perigee_SOURCE_DIR}/../src/PDNSolution_U_Mixed_3D.cpp
  ${perigee_SOURCE_DIR}/../src/PNonlinear_Seg_Solver.cpp
  ${perigee_SOURCE_DIR}/../src/PTime_Seg_Solver.cpp
  )

# 3. Postprocess cpp
SET( perigee_postprocess_lib_src
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_SOURCE_DIR}/../src/VisDataPrep_ALE_NS_3D.cpp
  ${perigee_SOURCE_DIR}/../src/VisDataPrep_Mixed_U_Hyperelastic_3D.cpp
  ${perigee_SOURCE_DIR}/../src/VisDataPrep_Mixed_FSI_3D.cpp
  ${perigee_SOURCE_DIR}/../src/VTK_Writer_Fluids_ALE_Tet4.cpp
  ${perigee_SOURCE_DIR}/../src/VTK_Writer_FSI_Tet4.cpp
  ${perigee_SOURCE_DIR}/../src/Post_error_tet4_ns.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  ${perigee_source}/Mesh/QuadPts_vis_tet4.cpp
  )

# ===============================================
# Finish setting my own source code
# ===============================================
# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../../include )
INCLUDE_DIRECTORIES( ../include )

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
TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( gmshIO gmsh_process.cpp )
ADD_EXECUTABLE( preprocess_fsi preprocess_fsi_tets.cpp)
ADD_EXECUTABLE( fsi_tet4_3d fsi_tet4_driver.cpp)
ADD_EXECUTABLE( prepostproc prepost_tets.cpp)

ADD_EXECUTABLE( vis_fsi vis_driver.cpp)
ADD_EXECUTABLE( vis_fsi_fluid vis_fsi_fluid.cpp)
ADD_EXECUTABLE( vis_fsi_solid vis_fsi_solid.cpp)
ADD_EXECUTABLE( vis_fsi_wss vis_fsi_wss.cpp)

# Link libraries
TARGET_LINK_LIBRARIES( gmshIO ${EXTRA_LINK_LIBS} perigee_preprocess )
TARGET_LINK_LIBRARIES( preprocess_fsi ${EXTRA_LINK_LIBS} perigee_preprocess )
TARGET_LINK_LIBRARIES( fsi_tet4_3d ${EXTRA_LINK_LIBS} perigee_analysis )
TARGET_LINK_LIBRARIES( prepostproc ${EXTRA_LINK_LIBS} perigee_preprocess )
TARGET_LINK_LIBRARIES( vis_fsi ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_fluid ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_solid ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_wss ${EXTRA_LINK_LIBS} perigee_postprocess )

#EOF
