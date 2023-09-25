// ==================================================================
// gmsh_process.cpp
//
// Code that handles the gmsh file.
//
// Date Created: July 1 2017
// Author: Ju Liu
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "fsi_2cube_hex_p2.msh";

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> update_quadratic_hex_IEN( 0 );

  GIO -> update_quadratic_hex_IEN( 1 );

  GIO -> write_quadratic_sur_vtu(0,0,true);
  GIO -> write_quadratic_sur_vtu(1,0,true);

  GIO -> write_quadratic_sur_vtu(2,0,true);
  GIO -> write_quadratic_sur_vtu(3,0,true);
  GIO -> write_quadratic_sur_vtu(4,0,true);
  GIO -> write_quadratic_sur_vtu(5,0,true);

  GIO -> write_quadratic_sur_vtu(5,1,true);
  GIO -> write_quadratic_sur_vtu(6,1,true);
  GIO -> write_quadratic_sur_vtu(7,1,true);
  GIO -> write_quadratic_sur_vtu(8,1,true);
  GIO -> write_quadratic_sur_vtu(9,1,true);
  GIO -> write_quadratic_sur_vtu(10,1,true);

  GIO -> write_each_vtu();

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO;
}

// EOF
