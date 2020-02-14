// ==================================================================
// gmsh_process.cpp
//
// Date Created: July 1 2017
// Author: Ju Liu
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "fsi_cylinder.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  // In FSI problems, we take domain tag 0 be the fluid domain;
  // domain tag 1 be the solid domain.
  GIO -> check_FSI_ordering();

  GIO -> print_info();

  // Write the fluid surface meshes associated with the fluid volume mesh
  GIO -> write_vtp(0,0,true);
  GIO -> write_vtp(1,0,true);
  GIO -> write_vtp(2,0);

  // Write the solid surface meshes associated with the solid volume mesh
  GIO -> write_vtp(3,1);
  GIO -> write_vtp(4,1);
  GIO -> write_vtp(5,1);

  // Write each 3D volumetric mesh separately
  GIO -> write_each_vtu();

  // Write the whole volumetric mesh with each sub-domain with a physical tag
  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
