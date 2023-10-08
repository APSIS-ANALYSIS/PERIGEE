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

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  const std::vector<int> nbc_face_id = config["nbc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> nbc_face_name = config["nbc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> nbc_vol_id = config["nbc_vol_id"].as<std::vector<int>>();

  const std::vector<int> ebc_face_id = config["ebc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> ebc_face_name = config["ebc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> ebc_vol_id = config["ebc_vol_id"].as<std::vector<int>>();

  const int eleType = GIO -> get_eleType(0);
  const int num_phy_domain_1d = GIO -> get_num_phy_domain_1d();
  const int num_phy_domain_2d = GIO -> get_num_phy_domain_2d();
  
  // check whether it is a 3d problem or not
  SYS_T::print_fatal_if( num_phy_domain_1d != 0 || num_phy_domain_2d <= 0,
    "ERROR: the number of 1D physical domain should be 0. \n" );

  if( eleType==2 || eleType==3 )
  {
    for( int ii=0; ii<nbc_face_id.size(); ++ii )
      GIO -> write_vtp( nbc_face_id[ii], nbc_vol_id[ii], false );

    for( int ii=0; ii<ebc_face_id.size(); ++ii )
      GIO -> write_vtp( ebc_face_id[ii], ebc_vol_id[ii], true );
  }
  else if( eleType==9 )
  {
    GIO -> update_quadratic_tet_IEN(0);
    for( int ii=0; ii<nbc_face_id.size(); ++ii )
      GIO -> write_quadratic_sur_vtu( nbc_face_id[ii], nbc_vol_id[ii], false );

    for( int ii=0; ii<ebc_face_id.size(); ++ii )
      GIO -> write_quadratic_sur_vtu( ebc_face_id[ii], ebc_vol_id[ii], true );
  }
  else if( eleType==10 )
  {
    GIO -> update_quadratic_hex_IEN(0);
    for( int ii=0; ii<nbc_face_id.size(); ++ii )
      GIO -> write_quadratic_sur_vtu( nbc_face_id[ii], nbc_vol_id[ii], false );

    for( int ii=0; ii<ebc_face_id.size(); ++ii )
      GIO -> write_quadratic_sur_vtu( ebc_face_id[ii], ebc_vol_id[ii], true );
  }
  else SYS_T::print_fatal("Error: the element type of gmsh file cannot be read. \n");

  const std::string wmname = config["wmname"].as<std::string>();
  const bool isXML = config["isXML"].as<bool>();
  GIO -> write_vtu( wmname, isXML );

  delete GIO;
  return EXIT_SUCCESS;
}
  
  

// EOF
