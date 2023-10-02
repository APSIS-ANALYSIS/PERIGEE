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
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::string gmshFile = "fsi_cylinder_hex_p1_wBL_test.msh";

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> write_vtp("lumen_inlet_vol_000",0,0,true);
  GIO -> write_vtp("lumen_outlet_vol_000",1,0,true);
  GIO -> write_vtp("lumen_wall_vol",2,0,true);

  GIO -> write_vtp("tissue_interior_wall_vol",2,1,true);
  GIO -> write_vtp("tissue_inlet_vol_000",3,1,true);
  GIO -> write_vtp("tissue_outlet_vol_000",4,1,true);
  GIO -> write_vtp("tissue_wall_vol",5,1,true);

  const std::vector<std::string> name_list {"lumen_vol","tissue_vol"};

  GIO -> write_each_vtu(name_list);

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}