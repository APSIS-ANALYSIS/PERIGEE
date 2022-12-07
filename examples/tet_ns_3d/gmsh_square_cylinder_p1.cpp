// ============================================================================
// gmsh_square_cylinder_p1.cpp
//
// This is the translator that transforms the square-cylinder geometry in P1
// element into vtk files.
//
// Date: Dec. 7 2022
// ============================================================================
#include "Gmsh_FileIO.hpp"

int main(int argc, char * argv[] )
{
  std::string gmshFile = "turbulent_cylinder_LR.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::GetOptionString("-gmsh_file", gmshFile);

  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
