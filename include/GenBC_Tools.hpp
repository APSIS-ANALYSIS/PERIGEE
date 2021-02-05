#ifndef GENBC_TOOLS_HPP
#define GENBC_TOOLS_hpp
// ==================================================================
// GenBC_Tools.hpp
//
// This is a suite of function tools that will be called in the reduced
// modeling of the cardiovascular system
// ==================================================================
#include "Sys_Tools.hpp"

namespace GENBC_T
{
  // ----------------------------------------------------------------
  // ! get_genbc_file_type : read the genbc file and determine what
  //   type of gen bc the file specifies. It will return
  //   0 for unknown type
  //   1 for Resistance
  //   2 for RCR
  //   3 for Inductance
  //   4 for Coronary (including RCR)
  // ----------------------------------------------------------------
  int get_genbc_file_type( const char * const &lpn_filename );
}

#endif
