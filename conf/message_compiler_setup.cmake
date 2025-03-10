# This file print the compiler information

# ===============================================
# Print Compiler/Linker and flags
# ===============================================
message(STATUS "=================== Compiler setup ==========================")
message(STATUS "C Compiler: \"${CMAKE_C_COMPILER}\"")
message(STATUS "C++ Compiler: \"${CMAKE_CXX_COMPILER}\"")
message(STATUS "C Compiler ID: \"${CMAKE_C_COMPILER_ID}\"")
message(STATUS "C++ Compiler ID: \"${CMAKE_CXX_COMPILER_ID}\"")
message(STATUS "C Flags: \"${CMAKE_C_FLAGS}\"")
message(STATUS "C++ Flags: \"${CMAKE_CXX_FLAGS}\"")
message(STATUS "Linker: \"${CMAKE_LINKER}\"")
message(STATUS "Executable Linker Flags: \"${CMAKE_EXE_LINKER_FLAGS}\"")
message(STATUS "Shared Linker Flags: \"${CMAKE_SHARED_LINKER_FLAGS}\"")

# ===============================================
# Print Misc Messages 
# ===============================================
message(STATUS "===================  Misc setup  ==========================")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
message(STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
message(STATUS "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS} )

message(STATUS "===========================================================")

# EOF
