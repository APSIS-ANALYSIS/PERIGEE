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
  ${perigee_source}/System/SV_Tools.cpp
  ${perigee_source}/Mesh/VTK_Tools.cpp
  )

# ===============================================
# Finish setting my own source code
# ===============================================
# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )
SET(YAML_DIR "${HOME_DIR}/lib/yaml-shared")
SET(YAML_INC "${YAML_DIR}/include/")
SET(YAML_LIB "${YAML_DIR}/lib/libyaml-cpp.so")
INCLUDE_DIRECTORIES(${YAML_INC})

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess ${EXTRA_LINK_LIBS} )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

FIND_PACKAGE(VTK COMPONENTS 
  CommonComputationalGeometry
  FiltersSources
  IOXML
)

IF (NOT VTK_FOUND)
  MESSAGE(FATAL_ERROR "gmsh_file_converter: Unable to find the VTK build folder.")
ENDIF()

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
ADD_EXECUTABLE( gmsh_file_converter MACOSX_BUNDLE gmsh_file_converter.cpp )

TARGET_LINK_LIBRARIES( gmsh_file_converter ${EXTRA_LINK_LIBS} perigee_preprocess ${VTK_LIBRARIES} ${YAML_LIB} )

# VTK_MODULE_AUTOINIT is needed
VTK_MODULE_AUTOINIT(
  TARGETS gmsh_file_converter
  MODULES ${VTK_LIBRARIES}
)

# EOF
