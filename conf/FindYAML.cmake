#############################################################
# Try to find YAML                                          #
#                                                           #
# Once done this will define:                               #
#  YAML_FOUND     - system has YAML                         #
#  YAML_DIR       - YAML directory                          #
#  YAML_INC       - YAML include directory                  #
#  YAML_LIB       - YAML library (static or dynamic)        #
#                                                           #
# Usage:                                                    #
#  find_package(YAML)                                       #
#                                                           #
# Setting these changes the behavior of the search          #
#  YAML_DIR       - YAML directory                          #
#############################################################

## Try to set YAML_DIR ##
##########################
if(NOT DEFINED YAML_DIR)
  set(YAML_DIR $ENV{YAML_DIR})
endif()

## Includes ##
##############
if(EXISTS "${YAML_DIR}/include")
 set(YAML_INC "${YAML_DIR}/include")
else()
  message(SEND_ERROR "YAML includes not found")
endif()

## Library ##
#############
if(EXISTS "${YAML_DIR}/lib/libyaml-cpp.so")
  set(YAML_LIB "${YAML_DIR}/lib/libyaml-cpp.so")
elseif(EXISTS "${YAML_DIR}/lib/libyaml-cpp.a")
  set(YAML_LIB "${YAML_DIR}/lib/libyaml-cpp.a")
elseif(EXISTS "${YAML_DIR}/lib/libyaml-cpp.dylib")
  set(YAML_LIB "${YAML_DIR}/lib/libyaml-cpp.dylib")
else()
  message(SEND_ERROR "YAML library not found")
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAML
  "YAML could not be found: be sure to set YAML_DIR in your environment variables"
  YAML_LIB YAML_INC YAML_DIR)

# EOF
