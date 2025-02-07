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

# print makefile details or not
SET(CMAKE_VERBOSE_MAKEFILE FALSE)

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)

# ===============================================
# Set my own source code
# ===============================================
SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
SET( perigee_lib_src 
  ${perigee_source}/System/Vector_3.cpp
  ${perigee_source}/System/SV_Tools.cpp
  ${perigee_source}/System/FE_Tools.cpp
  ${perigee_source}/Mesh/VTK_Tools.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  )

# ===============================================
# Finish setting my own source code
# ===============================================
# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )

# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee ${perigee_lib_src} )
TARGET_LINK_LIBRARIES( perigee ${EXTRA_LINK_LIBS} )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
# 1. Make executable
# 1.1 Preprocess
ADD_EXECUTABLE( sv_fsi_converter convert_sv_fsi_files.cpp)

# -------------------------------------------------------------------
# 2. Link libraries
TARGET_LINK_LIBRARIES( sv_fsi_converter ${EXTRA_LINK_LIBS} perigee )

#EOF
