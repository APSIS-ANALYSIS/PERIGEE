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
  std::string gmshFile = "fsi_2cube_hex_p1.msh";

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> write_vtp(0,0,true);
  GIO -> write_vtp(1,0,true);

  GIO -> write_vtp(2,0,true);
  GIO -> write_vtp(3,0,true);
  GIO -> write_vtp(4,0,true);
  GIO -> write_vtp(5,0,true);

  GIO -> write_vtp(5,1,true);
  GIO -> write_vtp(6,1,true);
  GIO -> write_vtp(7,1,true);
  GIO -> write_vtp(8,1,true);
  GIO -> write_vtp(9,1,true);
  GIO -> write_vtp(10,1,true);

  GIO -> write_each_vtu();

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
