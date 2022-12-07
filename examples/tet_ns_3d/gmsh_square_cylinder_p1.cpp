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

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();


  return EXIT_SUCCESS;
}


// EOF
