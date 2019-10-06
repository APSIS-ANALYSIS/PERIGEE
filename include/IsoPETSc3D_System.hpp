#ifndef ISOPETSC3D_SYSTEM_HPP
#define ISOPETSC3D_SYSTEM_HPP
// ==================================================================
// IsoPETSc3D_System.hpp
// The main IsoPETSc3D system include hpp file. 
// ==================================================================

#include <climits>
// control the int data type used int this project.
// 1. smart_int, could be switched to 64bit int when necessary
typedef int s_int;

// 2. little_int fixed to be 32bit int
typedef int l_int;

// 3. unsinged smart int: unsigned smart integer
typedef unsigned int us_int;

#endif
