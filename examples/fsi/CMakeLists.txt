CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(perigee)

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

# ===============================================
# Set my own source code
# ===============================================
SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/VTK_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/Hex_Tools.cpp
  ${perigee_source}/Mesh/ElemBC_3D.cpp
  ${perigee_source}/Mesh/ElemBC_3D_outflow.cpp
  ${perigee_source}/Mesh/NodalBC.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_1D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Quad.cpp
  ${perigee_source}/Mesh/NBC_Partition.cpp
  ${perigee_source}/Mesh/NBC_Partition_MF.cpp
  ${perigee_source}/Mesh/NBC_Partition_inflow.cpp
  ${perigee_source}/Mesh/NBC_Partition_inflow_MF.cpp
  ${perigee_source}/Mesh/EBC_Partition.cpp
  ${perigee_source}/Mesh/EBC_Partition_outflow_MF.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Mesh/Global_Part_Reload.cpp
  ${perigee_source}/Mesh/Part_FEM.cpp
  ${perigee_source}/Mesh/Part_FEM_FSI.cpp
  ${perigee_source}/Element/FE_Tools.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad4_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad9_3D_der0.cpp  
  ${perigee_SOURCE_DIR}/src/NodalBC_3D_FSI.cpp
  )

# 2. Anlysis source cpp
SET( perigee_analysis_lib_src
  ${perigee_source}/System/PETSc_Tools.cpp
  ${perigee_source}/System/Matrix_PETSc.cpp
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/Tensor2_3D.cpp
  ${perigee_source}/System/SymmTensor2_3D.cpp
  ${perigee_source}/System/Tensor4_3D.cpp
  ${perigee_source}/System/SymmTensor4_3D.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_1D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Quad.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Hex.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_NBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_InflowBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC_outflow.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/APart_Node_FSI.cpp
  ${perigee_source}/Analysis_Tool/FEANode.cpp
  ${perigee_source}/Element/FE_Tools.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  ${perigee_source}/Element/FEAElement_Tet10.cpp
  ${perigee_source}/Element/FEAElement_Hex8.cpp
  ${perigee_source}/Element/FEAElement_Hex27.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad4_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad9_3D_der0.cpp  
  ${perigee_source}/Model/Tissue_prestress.cpp
  ${perigee_source}/Model/GenBC_RCR.cpp
  ${perigee_source}/Model/GenBC_Resistance.cpp
  ${perigee_source}/Model/GenBC_Inductance.cpp
  ${perigee_source}/Model/GenBC_Coronary.cpp
  ${perigee_source}/Model/GenBC_Pressure.cpp
  ${perigee_source}/Model/GenBC_Tools.cpp
  ${perigee_source}/Model/FlowRate_Steady.cpp
  ${perigee_source}/Model/FlowRate_Unsteady.cpp
  ${perigee_source}/Model/FlowRate_Linear2Steady.cpp
  ${perigee_source}/Model/FlowRate_Cosine2Steady.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_FSI_Mesh_Elastostatic.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_FSI_Mesh_Laplacian.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_2x2Block_VMS_Incompressible.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_2x2Block_VMS_Hyperelasticity.cpp
  ${perigee_SOURCE_DIR}/src/PDNSolution_V.cpp
  ${perigee_SOURCE_DIR}/src/PDNSolution_P.cpp
  ${perigee_SOURCE_DIR}/src/PGAssem_Mesh.cpp
  ${perigee_SOURCE_DIR}/src/PGAssem_FSI.cpp
  # ${perigee_SOURCE_DIR}/src/PGAssem_Wall_Prestress.cpp
  ${perigee_SOURCE_DIR}/src/PNonlinear_FSI_Solver.cpp
  ${perigee_SOURCE_DIR}/src/PTime_FSI_Solver.cpp
  )

# 3. Postprocess cpp
SET( perigee_postprocess_lib_src
  ${perigee_source}/Mesh/VTK_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/Hex_Tools.cpp
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_FSI.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_ALE_NS.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_Hyperelastic.cpp
  ${perigee_SOURCE_DIR}/src/VTK_Writer_FSI.cpp
  )

# ===============================================
# Finish setting my own source code
# ===============================================
# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )
INCLUDE_DIRECTORIES( ./include )

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

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( preprocess_fsi preprocess.cpp)
# ADD_EXECUTABLE( wall_ps3d wall_ps_driver.cpp)
ADD_EXECUTABLE( fsi3d driver.cpp)
ADD_EXECUTABLE( prepostproc prepost.cpp)
ADD_EXECUTABLE( vis_fsi vis.cpp )
ADD_EXECUTABLE( vis_fluid vis_fluid.cpp )
ADD_EXECUTABLE( vis_solid vis_solid.cpp )
ADD_EXECUTABLE( vis_fsi_wss vis_fsi_wss.cpp )
ADD_EXECUTABLE( vis_fsi_wss_tet4 vis_fsi_wss_tet4.cpp )
ADD_EXECUTABLE( vis_fsi_wss_hex8 vis_fsi_wss_hex8.cpp )

# Link libraries
TARGET_LINK_LIBRARIES( preprocess_fsi ${EXTRA_LINK_LIBS} perigee_preprocess )
# TARGET_LINK_LIBRARIES( wall_ps3d ${EXTRA_LINK_LIBS} perigee_analysis )
TARGET_LINK_LIBRARIES( fsi3d ${EXTRA_LINK_LIBS} perigee_analysis )
TARGET_LINK_LIBRARIES( prepostproc ${EXTRA_LINK_LIBS} perigee_preprocess )
TARGET_LINK_LIBRARIES( vis_fsi ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fluid ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_solid ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_wss ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_wss_tet4 ${EXTRA_LINK_LIBS} perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_fsi_wss_hex8 ${EXTRA_LINK_LIBS} perigee_postprocess )

if(OPENMP_CXX_FOUND)
  SET_TARGET_PROPERTIES( perigee_preprocess PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  SET_TARGET_PROPERTIES( preprocess_fsi PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  SET_TARGET_PROPERTIES( prepostproc PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  TARGET_INCLUDE_DIRECTORIES( perigee_preprocess PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_INCLUDE_DIRECTORIES( preprocess_fsi PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_INCLUDE_DIRECTORIES( prepostproc PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_LINK_LIBRARIES( perigee_preprocess ${OpenMP_CXX_LIBRARIES} )
endif()

#EOF
