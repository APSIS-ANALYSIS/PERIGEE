#ifndef FLOWRATE_TOOLS_HPP
#define FLOWRATE_TOOLS_HPP

#include "Sys_Tools.hpp"

namespace FLOW_T
{
  // ----------------------------------------------------------------
  // ! get_flowrate_file_type : read the genbc file and determine what
  //   type of gen bc the file specifies. It will return
  //   0 : unknown type
  //   1 : Steady
  //   2 : Unsteady
  //   3 : Linear2Steady
  //   4 : Cosine2Steady
  // ----------------------------------------------------------------
  int get_flowrate_file_type( const std::string &inflow_filename );
}

#endif