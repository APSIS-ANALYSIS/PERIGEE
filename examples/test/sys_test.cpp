#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"

#include "Gmsh_FileIO.hpp"
#include "yaml-cpp/yaml.h"
#include "omp.h"

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
  // PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  // SYS_T::print_perigee_art();
  // SYS_T::print_system_info();
  // PetscFinalize();

  std::vector<std::string> YAML_Lists {"cube_tet_p1.yml", "cube_tet_p2.yml", "cube_hex_p1.yml", "cube_hex_p2.yml"};
  std::vector<std::string> ele_type_tag {"tet_p1_", "tet_p2_", "hex_p1_", "hex_p2_"};
  std::vector<std::vector<std::string>> FileLists {
    {"lef_vol.vtp", "rig_vol.vtp", "bac_vol.vtp", "top_vol.vtp", "bot_vol.vtp", "fro_vol.vtp", "whole_vol.vtu"},
    {"lef_vol.vtu", "rig_vol.vtu", "bac_vol.vtu", "top_vol.vtu", "bot_vol.vtu", "fro_vol.vtu", "whole_vol.vtu"},
    {"left_vol.vtp", "right_vol.vtp", "back_vol.vtp", "top_vol.vtp", "bottom_vol.vtp", "front_vol.vtp", "whole_vol.vtu"},
    {"left_vol.vtu", "right_vol.vtu", "back_vol.vtu", "top_vol.vtu", "bottom_vol.vtu", "front_vol.vtu", "whole_vol.vtu"}
  };

  for(int kk=0; kk<4; ++kk)
  {
    YAML::Node config = YAML::LoadFile(YAML_Lists[kk]);

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

    if( GIO->get_eleType(0)==2 || GIO->get_eleType(0)==3 )
    {
      for( int ii=0; ii<nbc_face_id.size(); ++ii )
        GIO -> write_vtp( nbc_face_id[ii], nbc_vol_id[ii], false );

      for( int ii=0; ii<ebc_face_id.size(); ++ii )
        GIO -> write_vtp( ebc_face_id[ii], ebc_vol_id[ii], true );
    }
    else if( GIO->get_eleType(0)==9 )
    {
      GIO -> update_quadratic_tet_IEN(0);
      for( int ii=0; ii<nbc_face_id.size(); ++ii )
        GIO -> write_quadratic_sur_vtu( nbc_face_id[ii], nbc_vol_id[ii], false );

      for( int ii=0; ii<ebc_face_id.size(); ++ii )
        GIO -> write_quadratic_sur_vtu( ebc_face_id[ii], ebc_vol_id[ii], true );
    }
    else if( GIO->get_eleType(0)==10 )
    {
      GIO -> update_quadratic_hex_IEN(0);
      for( int ii=0; ii<nbc_face_id.size(); ++ii )
        GIO -> write_quadratic_sur_vtu( nbc_face_id[ii], nbc_vol_id[ii], false );

      for( int ii=0; ii<ebc_face_id.size(); ++ii )
        GIO -> write_quadratic_sur_vtu( ebc_face_id[ii], ebc_vol_id[ii], true );
    }
    else
      { SYS_T::print_fatal("Error: the element type of gmsh file cannot be read. \n"); return 0; }

    const std::string wmname = config["wmname"].as<std::string>();
    const bool isXML = config["isXML"].as<bool>();
    GIO -> write_vtu( wmname, isXML );

    delete GIO;
    for(int jj=0; jj<7; ++jj)
    {
      std::string temp = ele_type_tag[kk];
      if( compare(FileLists[kk][jj], temp.append(FileLists[kk][jj])) )
        std::cout<<ele_type_tag[kk]+FileLists[kk][jj]<<" are IDENTICAL!"<<std::endl;
      else
        std::cout<<ele_type_tag[kk]+FileLists[kk][jj]<<" are NOT IDENTICAL!"<<std::endl;
    }
    std::cout<<std::endl;
  }
  
  #pragma omp parallel
  {
    SYS_T::print_fatal("TEST PRINT_FATAL WITH OPENMP.");
  }

  return EXIT_SUCCESS;
}

// EOF
