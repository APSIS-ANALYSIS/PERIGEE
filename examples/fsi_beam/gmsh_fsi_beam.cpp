// ==================================================================
// gmsh_beam_process.cpp
//
// Code that handles the gmsh file for a patient-specific aorta 
// geometry.
//
// Date Created: Oct. 12  2017
// Author: Ju Liu
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "fsi_beam.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  // In FSI problems, we take domain tag 0 be the fluid domain;
  // domain tag 1 be the solid domain.
  GIO -> check_FSI_ordering();

  GIO -> print_info();

  // Write the fluid surface meshes associated with the fluid volume mesh
  // Inlet face
  // I need to set this as true to associate with the volume mesh, since
  // I need this in the calculation of the outward normal vector
  
  // 0 is fluid; 1 is solid
  GIO -> write_vtp(0,0);
  GIO -> write_vtp(1,0);
  GIO -> write_vtp(2,0);
  GIO -> write_vtp(3,0);
  GIO -> write_vtp(4,0);
  GIO -> write_vtp(5,0, true);
  GIO -> write_vtp(6,0);
  
  // Write the solid surface meshes associated with the solid volume mesh
  GIO -> write_vtp(7,1); // solid wall
  GIO -> write_vtp(8,1); // solid caps
  GIO -> write_vtp(9,1); // solid caps
  GIO -> write_vtp(10,1); // solid caps
  GIO -> write_vtp(11,1); // solid caps
  GIO -> write_vtp(12,1); // solid caps

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
