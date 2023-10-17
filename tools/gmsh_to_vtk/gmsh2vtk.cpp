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
    // Set number of threads
#ifdef _OPENMP
  omp_set_num_threads( omp_get_num_procs() );
  SYS_T::print_omp_info();
#endif

  const std::string input_yaml_file("example.yml");
  
  SYS_T::print_fatal_if( !SYS_T::file_exist( input_yaml_file ), "ERROR: the file %s does not exist on disk.\n", input_yaml_file.c_str() );

  // Get element type from GIO
  YAML::Node config = YAML::LoadFile( input_yaml_file );

  const std::string gmshFile = config["gmsh_file"].as<std::string>();
  const double isFSI = config["isFSI"].as<bool>();

  SYS_T::print_fatal_if( !SYS_T::file_exist( gmshFile ), "ERROR: the file %s does not exist on disk.\n", gmshFile.c_str() );

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  if( isFSI ) GIO -> check_FSI_ordering();

  GIO -> print_info();

  const std::vector<int> nbc_face_id = config["nbc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> nbc_face_name = config["nbc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> nbc_vol_id = config["nbc_vol_id"].as<std::vector<int>>();
  const std::vector<std::string> nbc_face_file_name = config["nbc_face_file_name"].as<std::vector<std::string>>();
  const std::vector<bool> nbc_isXML = config["nbc_isXML"].as<std::vector<bool>>();

  const std::vector<int> ebc_face_id = config["ebc_face_id"].as<std::vector<int>>();
  const std::vector<std::string> ebc_face_name = config["ebc_face_name"].as<std::vector<std::string>>();
  const std::vector<int> ebc_vol_id = config["ebc_vol_id"].as<std::vector<int>>();
  const std::vector<std::string> ebc_face_file_name = config["ebc_face_file_name"].as<std::vector<std::string>>();
  const std::vector<bool> ebc_isXML = config["ebc_isXML"].as<std::vector<bool>>();

  const int eleType = GIO -> get_eleType(0);
  const int num_phy_domain_1d = GIO -> get_num_phy_domain_1d();
  const int num_phy_domain_2d = GIO -> get_num_phy_domain_2d();
  
  // check whether it is a 3d problem or not
  SYS_T::print_fatal_if( num_phy_domain_1d != 0 || num_phy_domain_2d <= 0,
    "ERROR: the number of 1D physical domain should be 0. \n" );

  if( eleType==2 || eleType==3 )
  {
    for( unsigned int ii=0; ii<nbc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== nbc_face_name: "<<nbc_face_name[ii]<<'\t'; 
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(nbc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(nbc_vol_id[ii])<<std::endl;
       
      GIO -> write_vtp( nbc_face_file_name[ii], nbc_face_id[ii], nbc_vol_id[ii], nbc_isXML[ii] );
    }

    for( unsigned int ii=0; ii<ebc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== ebc_face_name: "<<ebc_face_name[ii]<<'\t';
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(ebc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(ebc_vol_id[ii])<<std::endl;
  
      GIO -> write_vtp( ebc_face_file_name[ii], ebc_face_id[ii], ebc_vol_id[ii], ebc_isXML[ii] );
    }
  }
  else if( eleType==9 )
  {
    GIO -> update_quadratic_tet_IEN(0);
    if( isFSI ) GIO -> update_quadratic_tet_IEN(1);
    for( unsigned int ii=0; ii<nbc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== nbc_face_name: "<<nbc_face_name[ii]<<'\t'; 
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(nbc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(nbc_vol_id[ii])<<std::endl;

      GIO -> write_quadratic_sur_vtu( nbc_face_file_name[ii], nbc_face_id[ii], nbc_vol_id[ii], nbc_isXML[ii] );
    }
    for( unsigned int ii=0; ii<ebc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== ebc_face_name: "<<ebc_face_name[ii]<<'\t';
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(ebc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(ebc_vol_id[ii])<<std::endl;

      GIO -> write_quadratic_sur_vtu( ebc_face_file_name[ii], ebc_face_id[ii], ebc_vol_id[ii], ebc_isXML[ii] );
    }
  }
  else if( eleType==10 )
  {
    GIO -> update_quadratic_hex_IEN(0);
    if( isFSI ) GIO -> update_quadratic_hex_IEN(1);
    for( unsigned int ii=0; ii<nbc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== nbc_face_name: "<<nbc_face_name[ii]<<'\t'; 
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(nbc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(nbc_vol_id[ii])<<std::endl;

      GIO -> write_quadratic_sur_vtu( nbc_face_file_name[ii], nbc_face_id[ii], nbc_vol_id[ii], nbc_isXML[ii] );
    }
    for( unsigned int ii=0; ii<ebc_face_id.size(); ++ii )
    {
      std::cout<<'\n'<<"=== ebc_face_name: "<<ebc_face_name[ii]<<'\t';
      std::cout<<"phy_face_name: "<<GIO->get_phy_name_2d(ebc_face_id[ii])<<" with ";
      std::cout<<GIO->get_phy_name_3d(ebc_vol_id[ii])<<std::endl;

      GIO -> write_quadratic_sur_vtu( ebc_face_file_name[ii], ebc_face_id[ii], ebc_vol_id[ii], ebc_isXML[ii] );
    }  
  }
  else SYS_T::print_fatal("Error: the element type of gmsh file cannot be read. \n");

  if( isFSI ) GIO -> write_each_vtu();
  const std::string wmname = config["wmname"].as<std::string>();
  const bool isXML = config["isXML"].as<bool>();
  GIO -> write_vtu( wmname, isXML );

  delete GIO;
  return EXIT_SUCCESS;
}

// EOF
