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
  std::string gmshFile = "cylinder_p2.msh";
  int num_outlet = 1;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;
  std::cout<<" -num_outlet: "<<num_outlet<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  GIO -> update_quadratic_tet_IEN( 0 );

  GIO -> write_quadratic_sur_vtu( 0, 0, false); // assumed to be wall vtp
  GIO -> write_quadratic_sur_vtu( 1, 0, true); // assumed to be wall vtp
  
  for(int ii=0; ii<num_outlet; ++ii)
    GIO -> write_quadratic_sur_vtu( 2+ii, 0, true );

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
