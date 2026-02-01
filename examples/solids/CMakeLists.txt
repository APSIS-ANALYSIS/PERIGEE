CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

project(perigee)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

SET(perigee_VERSION_MAJOR 0)
SET(perigee_VERSION_MINOR 0)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "set to release" FORCE)
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "default to debug" FORCE)
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

OPTION(BUILD_SHARED_LIBS "build shared libraries" OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)

SET( perigee_source ${perigee_SOURCE_DIR}/../../src )

SET( PERIGEE_INCLUDE_DIRS 
  ${perigee_SOURCE_DIR}/include 
  ${perigee_SOURCE_DIR}/../../include )

# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/VTK_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/Hex_Tools.cpp
  ${perigee_source}/Mesh/Part_FEM.cpp
  ${perigee_source}/Mesh/ElemBC_3D.cpp
  ${perigee_source}/Mesh/ElemBC_3D_outflow.cpp
  ${perigee_source}/Mesh/ElemBC_3D_WallModel.cpp
  ${perigee_source}/Mesh/NodalBC.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_1D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Quad.cpp
  ${perigee_source}/Mesh/NBC_Partition.cpp
  ${perigee_source}/Mesh/NBC_Partition_inflow.cpp
  ${perigee_source}/Mesh/EBC_Partition.cpp
  ${perigee_source}/Mesh/EBC_Partition_outflow.cpp
  ${perigee_source}/Mesh/EBC_Partition_WallModel.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Element/FE_Tools.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad4_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Quad9_3D_der0.cpp  
  )

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
  ${perigee_source}/Mesh/QuadPts_Gauss_1D.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Hex.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Quad.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_NBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_InflowBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC_outflow.cpp
  ${perigee_source}/Analysis_Tool/ALocal_WeakBC.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
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
  ${perigee_source}/Model/GenBC_RCR.cpp
  ${perigee_source}/Model/GenBC_Resistance.cpp
  ${perigee_source}/Model/GenBC_Inductance.cpp
  ${perigee_source}/Model/GenBC_Coronary.cpp
  ${perigee_source}/Model/GenBC_Tools.cpp
  ${perigee_source}/Model/GenBC_Pressure.cpp
  ${perigee_source}/Model/FlowRate_Steady.cpp
  ${perigee_source}/Model/FlowRate_Unsteady.cpp
  ${perigee_source}/Model/FlowRate_Linear2Steady.cpp
  ${perigee_source}/Model/FlowRate_Cosine2Steady.cpp
  ${perigee_source}/Model/FlowRate_Sine2Zero.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/src/PNonlinear_Solid_Solver.cpp
  ${perigee_SOURCE_DIR}/src/PTime_Solid_Solver.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_2x2Block_VMS_Hyperelasticity.cpp
  ${perigee_SOURCE_DIR}/src/PGAssem_Solid_FEM.cpp
  )

SET( perigee_postprocess_lib_src
  ${perigee_source}/Mesh/VTK_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/Mesh/Hex_Tools.cpp
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  )

INCLUDE_DIRECTORIES( ${PERIGEE_INCLUDE_DIRS} )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess PUBLIC ${EXTRA_LINK_LIBS} )

# 2. Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis PUBLIC ${EXTRA_LINK_LIBS} )

# 3. Postprocess lib
ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess PUBLIC perigee_analysis )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( preprocess3d preprocess.cpp)
ADD_EXECUTABLE( solid3d driver.cpp)
#ADD_EXECUTABLE( ns3dherk ns_herk_driver.cpp)
#ADD_EXECUTABLE( ns3dherkA ns_herk_driver_accurateA.cpp)
ADD_EXECUTABLE( prepost3d prepost.cpp)
#ADD_EXECUTABLE( vis_ns vis.cpp)
#ADD_EXECUTABLE( vis_wss_tet4 vis_wss_tet4.cpp)
#ADD_EXECUTABLE( vis_wss_tet10 vis_wss_tet10.cpp)
#ADD_EXECUTABLE( vis_wss_hex8 vis_wss_hex8.cpp)
#ADD_EXECUTABLE( vis_wss_hex27 vis_wss_hex27.cpp)

TARGET_LINK_LIBRARIES( preprocess3d perigee_preprocess )
TARGET_LINK_LIBRARIES( solid3d perigee_analysis )
#TARGET_LINK_LIBRARIES( ns3dherk perigee_analysis )
#TARGET_LINK_LIBRARIES( ns3dherkA perigee_analysis )
TARGET_LINK_LIBRARIES( prepost3d perigee_preprocess )
#TARGET_LINK_LIBRARIES( vis_ns perigee_postprocess )
#TARGET_LINK_LIBRARIES( vis_wss_tet4 perigee_postprocess )
#TARGET_LINK_LIBRARIES( vis_wss_tet10 perigee_postprocess )
#TARGET_LINK_LIBRARIES( vis_wss_hex8 perigee_postprocess )
#TARGET_LINK_LIBRARIES( vis_wss_hex27 perigee_postprocess )


if(OPENMP_CXX_FOUND)
  SET_TARGET_PROPERTIES( perigee_preprocess PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  SET_TARGET_PROPERTIES( preprocess3d PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  SET_TARGET_PROPERTIES( prepost3d PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  TARGET_INCLUDE_DIRECTORIES( perigee_preprocess PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_INCLUDE_DIRECTORIES( preprocess3d PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_INCLUDE_DIRECTORIES( prepost3d PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_LINK_LIBRARIES( perigee_preprocess PRIVATE ${OpenMP_CXX_LIBRARIES} )
endif()

# EOF
