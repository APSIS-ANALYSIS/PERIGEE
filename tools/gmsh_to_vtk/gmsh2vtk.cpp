// ============================================================================
// gmsh2vtk.cpp
//
// Code that handles the gmsh msh file and write the volumetric and surface data
// into vtk files.
// ============================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "cylinder.msh";
  int num_outlet = 1;

  for(int ii=1; ii<argc; ++ii)
  {
    if( std::string(argv[ii])=="-gmsh_file"  && ii+1<argc) gmshFile=argv[ii+1];
    if( std::string(argv[ii])=="-num_outlet" && ii+1<argc) num_outlet=std::stoi(argv[ii+1]);   
  }

  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;
  std::cout<<" -num_outlet: "<<num_outlet<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  delete GIO; 
  return EXIT_SUCCESS;
}

// EOF
