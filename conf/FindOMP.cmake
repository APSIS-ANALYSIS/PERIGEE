#############################################################
# Try to find OMP                                           #
#                                                           #
# Once done this will define:                               #
#   OMP_FOUND      - system has OMP                         #
#   OMP_DIR        - OMP directory                          #
#   OMP_INC        - OMP include directory                  #
#   OMP_LIB        - OMP library (static or dynamic)        #
#                                                           #
# Usage:                                                    #
#  find_package(OMP)                                        #
#                                                           #
# Setting these changes the behavior of the search          #
#  OMP_DIR         - OMP directory                          #
#############################################################

## Try to set OMP_DIR   ##
##########################
if(NOT DEFINED OMP_DIR)
	set(OMP_DIR $ENV{OMP_DIR})
endif()

## Includes ##
##############
if(EXISTS "${OMP_DIR}/include")
	set(OMP_INC "${OMP_DIR}/include")
else()
	message(SEND_ERROR "OpenMP includes not found")
endif()

## Library ##
#############
if(EXISTS "${OMP_DIR}/lib/libomp.so")
	set(OMP_LIB "${OMP_DIR}/lib/libomp.so")
elseif(EXISTS "${OMP_DIR}/lib/libomp.a")
	set(OMP_LIB "${OMP_DIR}/lib/libomp.a")
elseif(EXISTS "${OMP_DIR}/lib/libomp.dylib")
	set(OMP_LIB "${OMP_DIR}/lib/libomp.dylib")
else()
	message(SEND_ERROR "OpenMP library not found")
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OMP
	"OpenMP could not be found: be sure to set OMP_DIR in your environment variables"
	OMP_LIB OMP_INC OMP_DIR)

# EOF
