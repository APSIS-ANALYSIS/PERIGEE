// ============================================================================
// gmsh2vtk.cpp
//
// Code that handles the gmsh msh file and write the volumetric and surface data
// into vtk files.
// ============================================================================
#include "Gmsh_FileIO.hpp"
#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Get element type from GIO

  YAML::Node config = YAML::LoadFile("example.yml");

  const std::string gmshFile = config["gmsh_file"].as<std::string>();
  const int num_outlet = config["num_outlet"].as<int>();

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  const std::vector<int> nbc_face_id = config["nbc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> nbc_face_name = config["nbc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> nbc_vol_id = config["nbc_vol_id"].as<std::vector<int>>();

  const std::vector<int> ebc_face_id = config["ebc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> ebc_face_name = config["ebc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> ebc_vol_id = config["ebc_vol_id"].as<std::vector<int>>();

  // if element type is linear tet or hex
  // NEED TO ADD FILENAME TO BE WRITTEN
  for( int ii=0; ii<nbc_face_id.size(); ++ii )
    GIO -> write_vtp( nbc_face_id[ii], nbc_vol_id[ii], false );

  for( int ii=0; ii<ebc_face_id.size(); ++ii )
    GIO -> write_vtp( ebc_face_id[ii], ebc_vol_id[ii], true );

  // else if element type is quadratic tet or quadratic hex
  // ... update node numbering and write into vtu (not vtp)
  // else print message

  const std::string wmname = config["wmname"].as<std::string>();
  const bool isXML = config["isXML"].as<bool>();
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  return EXIT_SUCCESS;
}

// EOF
