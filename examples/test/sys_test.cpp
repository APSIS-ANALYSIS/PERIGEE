#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"

#include "Gmsh_FileIO.hpp"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>
#include <string>

bool compare(std::string& file1Path, std::string& file2Path) 
{
  std::ifstream file1(file1Path);
  std::ifstream file2(file2Path);

  if (!file1.is_open() || !file2.is_open())
  {
    std::cerr << "Failed to open one or both files." << std::endl;
    return false;
  }

  std::string line1, line2;
  int lineNumber = 1;
  while (std::getline(file1, line1) && std::getline(file2, line2))
  {
    if (line1 != line2) return false;
    lineNumber++;
  }

  // Check if one file is longer than the other
  if (std::getline(file1, line1) || std::getline(file2, line2)) return false;

  return true;
}


int main(int argc, char *argv[])
{
  std::vector<std::string> YAML_Lists {"cube_tet_p1.yml", "cube_tet_p2.yml", "cube_hex_p1.yml", "cube_hex_p2.yml"};
  std::vector<std::string> ele_type_tag {"tet_p1_", "tet_p2_", "hex_p1_", "hex_p2_"}; 
  std::vector<std::vector<std::string>> FileLists {
      {"lef.vtp", "rig.vtp", "bac.vtp", "top.vtp", "bot.vtp", "fro.vtp", "whole.vtu"},
      {"lef.vtu", "rig.vtu", "bac.vtu", "top.vtu", "bot.vtu", "fro.vtu", "whole.vtu"},
      {"lef.vtp", "rig.vtp", "bac.vtp", "top.vtp", "bot.vtp", "fro.vtp", "whole.vtu"},
      {"lef.vtu", "rig.vtu", "bac.vtu", "top.vtu", "bot.vtu", "fro.vtu", "whole.vtu"}
  };

  for(int kk=0; kk<4; ++kk)
  {
    const std::string input_yaml_file( YAML_Lists[kk] );
  
    SYS_T::print_fatal_if( !SYS_T::file_exist( input_yaml_file ), 
      "ERROR: the file %s does not exist on disk.\n", input_yaml_file.c_str() );

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
    for(int jj=0; jj<7; ++jj)
    {
      std::string temp = ele_type_tag[kk];
      if( compare(FileLists[kk][jj], temp.append(FileLists[kk][jj]) ) )
        std::cout<<ele_type_tag[kk]+FileLists[kk][jj]<<" are IDENTICAL!"<<std::endl;
      else
        std::cout<<ele_type_tag[kk]+FileLists[kk][jj]<<" are NOT IDENTICAL!"<<std::endl;
    }
    std::cout<<std::endl;
  }

  PERIGEE_OMP_PARALLEL
  {
    SYS_T::print_fatal("TEST PRINT_FATAL WITH OPENMP.");
  }

  return EXIT_SUCCESS;
}

// EOF
