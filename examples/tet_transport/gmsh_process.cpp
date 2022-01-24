// ============================================================================
// gmsh_process.cpp
//
// Handle the .msh file.
// ============================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "cylinder.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  const int vol_id = 0; 
  GIO -> write_vtp( 0, vol_id, true); 
  GIO -> write_vtp( 1, vol_id, true); 
  GIO -> write_vtp( 2, vol_id, true); 
  GIO -> write_vtp( 3, vol_id, true); 
  GIO -> write_vtp( 4, vol_id, true); 
  GIO -> write_vtp( 5, vol_id, true); 
  
  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
