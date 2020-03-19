# This file print the compiler information

# ===============================================
# Print Compiler/Linker and flags
# ===============================================
message(STATUS "=================== Compiler setup ==========================")
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})
message(STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER})
message(STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "CMAKE_COMPILER_IS_GNUCC: " ${CMAKE_COMPILER_IS_GNUCC})
message(STATUS "CMAKE_COMPILER_IS_GNUCXX: " ${CMAKE_COMPILER_IS_GNUCXX})

# ===============================================
# Print Misc Messages 
# ===============================================
message(STATUS "=================== Misc setup ==========================")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
message(STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
message(STATUS "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS} )

# EOF
