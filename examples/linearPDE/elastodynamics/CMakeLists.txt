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

INCLUDE(../../../conf/basic_variable_log.cmake)
INCLUDE(../../../conf/system_lib_loading.cmake)

SET( perigee_source ${perigee_SOURCE_DIR}/../../../src )

SET( PERIGEE_INCLUDE_DIRS 
  ${perigee_SOURCE_DIR}/../include 
  ${perigee_SOURCE_DIR}/../../../include )

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
  ${perigee_source}/Mesh/NodalBC.cpp
  ${perigee_source}/Mesh/NBC_Partition.cpp
  ${perigee_source}/Mesh/EBC_Partition.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Element/FE_Tools.cpp  
)

SET( perigee_analysis_lib_src
  ${perigee_source}/System/PETSc_Tools.cpp
  ${perigee_source}/System/Matrix_PETSc.cpp
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/Tensor2_3D.cpp
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
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/../src/PDNSolution_Elastodynamics.cpp
  ${perigee_SOURCE_DIR}/../src/PLocAssem_Elastodynamics_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/../src/PGAssem_LinearPDE_GenAlpha.cpp
  ${perigee_SOURCE_DIR}/../src/PNonlinear_LinearPDE_Solver.cpp
  ${perigee_SOURCE_DIR}/../src/PTime_LinearPDE_Solver.cpp
  )

SET( perigee_postprocess_lib_src
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  ${perigee_SOURCE_DIR}/../src/Post_error_elastodynamics.cpp
  ${perigee_SOURCE_DIR}/../src/VisDataPrep_Elastodynamics.cpp
  ${perigee_SOURCE_DIR}/../src/VTK_Writer_Elastodynamics.cpp
  )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess PUBLIC
  ${YAML_LIB} ${PETSC_LIBRARIES} ${VTK_LIBRARIES}
  ${HDF5_LIBRARIES} ${METIS_LIBRARIES} )
TARGET_INCLUDE_DIRECTORIES( perigee_preprocess PUBLIC
  ${perigee_SOURCE_DIR}/../../../include ${HDF5_INCLUDE_DIRS}
  ${PETSC_INCLUDE_DIRS} ${METIS_INCLUDE_DIR} )

# 2. Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis PUBLIC
  ${PETSC_STATIC_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARIES} )
TARGET_INCLUDE_DIRECTORIES( perigee_analysis PUBLIC
  ${PERIGEE_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}
  ${PETSC_INCLUDE_DIRS} )

# 3. Postprocess lib
ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess PUBLIC perigee_analysis )
TARGET_INCLUDE_DIRECTORIES( perigee_postprocess PUBLIC ${PERIGEE_INCLUDE_DIRS} )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( preprocess3d ../preprocess.cpp)
ADD_EXECUTABLE( elasto_3d driver.cpp)
ADD_EXECUTABLE( prepost3d ../prepost.cpp)
ADD_EXECUTABLE( post_error post_compare_manu.cpp )
ADD_EXECUTABLE( vis_elasto_3d vis.cpp )

TARGET_INCLUDE_DIRECTORIES( preprocess3d PRIVATE ${YAML_INC} )
TARGET_INCLUDE_DIRECTORIES( prepost3d PRIVATE ${YAML_INC} )

TARGET_LINK_LIBRARIES( preprocess3d PRIVATE perigee_preprocess )
TARGET_LINK_LIBRARIES( elasto_3d PRIVATE perigee_analysis )
TARGET_LINK_LIBRARIES( prepost3d PRIVATE perigee_preprocess )
TARGET_LINK_LIBRARIES( post_error PRIVATE perigee_postprocess )
TARGET_LINK_LIBRARIES( vis_elasto_3d PRIVATE perigee_postprocess )

if(OPENMP_CXX_FOUND)
  SET_TARGET_PROPERTIES( perigee_preprocess PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  SET_TARGET_PROPERTIES( preprocess3d PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  TARGET_INCLUDE_DIRECTORIES( perigee_preprocess PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_INCLUDE_DIRECTORIES( preprocess3d PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  TARGET_LINK_LIBRARIES( perigee_preprocess ${OpenMP_CXX_LIBRARIES} )
  TARGET_LINK_LIBRARIES( preprocess3d ${OpenMP_CXX_LIBRARIES} )
endif()

# EOF
