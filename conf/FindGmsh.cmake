#############################################################
# Try to find Gmsh                                          #
#                                                           #
# Once done this will define:                               #
#  GMSH_FOUND     - system has Gmsh                         #
#  GMSH_DIR       - Gmsh directory                          #
#  GMSH_INC       - Gmsh include directory                  #
#  GMSH_LIB       - Gmsh library (static or dynamic)        #
#                                                           #
# Usage:                                                    #
#  find_package( GMSH )                                     #
#                                                           #
# Setting these changes the behavior of the search          #
#  GMSH_DIR       - Gmsh directory                        #
#############################################################

## Try to set GMSH_DIR ##
##########################
if(NOT DEFINED GMSH_DIR)
  set(GMSH_DIR $ENV{GMSH_DIR})
endif()

## Includes ##
##############
if(EXISTS "${GMSH_DIR}/include")
 set(GMSH_INC "${GMSH_DIR}/include")
else()
  message(SEND_ERROR "GMSH includes not found")
endif()

## Library ##
#############
if(EXISTS "${GMSH_DIR}/lib/libgmsh.so")
  set(GMSH_LIB "${GMSH_DIR}/lib/libgmsh.so")
elseif(EXISTS "${GMSH_DIR}/lib/libgmsh.a")
  set(GMSH_LIB "${GMSH_DIR}/lib/libgmsh.a")
elseif(EXISTS "${GMSH_DIR}/lib/libgmsh.dylib")
  set(GMSH_LIB "${GMSH_DIR}/lib/libgmsh.dylib")
else()
  message(SEND_ERROR "GMSH library not found")
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  "GMSH could not be found: be sure to set GMSH_DIR in your environment variables"
  GMSH_LIB GMSH_INC GMSH_DIR)
