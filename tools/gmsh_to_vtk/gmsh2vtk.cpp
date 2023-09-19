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

  std::vector<int> nbc_face_id { 1, 2 };
  std::vector<std::string> nbc_face_name { "ftop", "fbot" };
  std::vector<int> nbc_vol_id { 0, 0 };

  std::vector<int> ebc_face_id { 0, 3, 4, 5 };
  std::vector<std::string> ebc_face_name { "ffro", "fbac", "flef", "frig" };
  std::vector<int> ebc_vol_id { 0, 0, 0, 0 };
  
  const int vol_id = 0;

  for( auto id : nbc_face_id )
    GIO -> write_vtp( id, vol_id, false );

  for( auto id : ebc_face_id )
    GIO -> write_vtp( id, vol_id, true );

  const std::string wmname("whole_vol");
  const bool ixXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  return EXIT_SUCCESS;
}

// EOF
