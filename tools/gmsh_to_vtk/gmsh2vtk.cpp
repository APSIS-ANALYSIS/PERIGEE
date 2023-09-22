// ============================================================================
// gmsh2vtk.cpp
//
// Code that handles the gmsh msh file and write the volumetric and surface data
// into vtk files.
// ============================================================================
#include "Gmsh_FileIO.hpp"
#include <yaml-cpp/yaml.h>

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

  YAML::Node config = YAML::LoadFile("example.yml");

  std::vector<int> nbc_face_id = config["nbc_face_id"].as<std::vector<int>>();
  std::vector<std::string> nbc_face_name = config["nbc_face_name"].as<std::vector<std::string>>();
  std::vector<int> nbc_vol_id = config["nbc_vol_id"].as<std::vector<int>>();

  std::vector<int> ebc_face_id = config["ebc_face_id"].as<std::vector<int>>();
  std::vector<std::string> ebc_face_name = config["ebc_face_name"].as<std::vector<std::string>>();
  std::vector<int> ebc_vol_id = config["ebc_vol_id"].as<std::vector<int>>();

  for( int ii=0; ii<nbc_face_id.size(); ++ii )
    GIO -> write_vtp( nbc_face_id[ii], nbc_vol_id[ii], false );

  for( int ii=0; ii<ebc_face_id.size(); ++ii )
    GIO -> write_vtp( ebc_face_id[ii], ebc_vol_id[ii], true );

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  return EXIT_SUCCESS;
}

// EOF
