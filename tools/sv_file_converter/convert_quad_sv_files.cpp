// ==================================================================
// This is a file convert function that correct the quadratic tet SV 
// files node and element indices to make them start from zero. The 
// default starting index in SimVascular is 1.
//
// Author: Ju Liu
// Date: Jan 20 2020
// ==================================================================
#include "SV_Tools.hpp"

int main( int argc, char * argv[] )
{
  // We assume the files contain a vtu file for the volume mesh
  // a single inlet and a single exterior wall mesh in vtp format
  // and several outlet mesh in vtp format, named as: outlet_xxx.vtp
  // Users are allowed to input the geo_file, sur_file_in, and
  // sur_file_wall, as well as sur_file_out_base. 
  std::string geo_file("vol.vtu");
  std::string sur_file_in("inlet.vtp");
  std::string sur_file_wall("wall.vtp");
  std::string sur_file_out_base("outlet_");
  int num_outlet = 1;

  std::string geo_out_name("whole_vol.vtu");
  std::string inl_out_name("inflow_vol.vtp");
  std::string wal_out_name("wall_vol.vtp");
  std::string out_out_base("outflow_vol_");

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  // Make sure this is a serial run
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1,"ERROR: converter is a serial routine! \n");
  
  SYS_T::GetOptionString("-vol_mesh", geo_file);
  SYS_T::GetOptionString("-inlet_mesh", sur_file_in);
  SYS_T::GetOptionString("-wall_mesh", sur_file_wall);
  SYS_T::GetOptionString("-outlet_mesh", sur_file_out_base);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  
  SYS_T::GetOptionString("-vol_geo_name", geo_out_name);
  SYS_T::GetOptionString("-inl_geo_name", inl_out_name);
  SYS_T::GetOptionString("-wal_geo_name", wal_out_name);
  SYS_T::GetOptionString("-out_geo_base", out_out_base);
  
  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -vol_mesh:     "<<geo_file<<std::endl; 
  std::cout<<" -inlet_mesh:   "<<sur_file_in<<std::endl; 
  std::cout<<" -wall_mesh:    "<<sur_file_wall<<std::endl; 
  std::cout<<" -outlet_mesh:  "<<sur_file_out_base<<std::endl; 
  std::cout<<" -num_outlet:   "<<num_outlet<<std::endl;
  std::cout<<"names to be written: \n";
  std::cout<<" -vol_geo_name: "<<geo_out_name<<std::endl;
  std::cout<<" -inl_geo_name: "<<inl_out_name<<std::endl;
  std::cout<<" -wal_geo_name: "<<wal_out_name<<std::endl;
  std::cout<<" -out_geo_base: "<<out_out_base<<std::endl;
  std::cout<<"================================="<<std::endl;

  // Boundary check for the number of outlets
  SYS_T::print_fatal_if(num_outlet < 1, 
      "Error: -num_outlet cannot be less than 1.\n");

  SYS_T::print_fatal_if(num_outlet > 1000, 
      "Error: -num_outlet cannot be more than 1000.\n");

  // If the corresponding name are the same, old file will be overwirtten
  // make a check.
  if( geo_file == geo_out_name || sur_file_in == inl_out_name ||
   sur_file_wall == wal_out_name || sur_file_out_base == out_out_base )
    std::cout<<"Warning: Some files have the same read and write names, new file will replace the old files.\n";

  // Generate the outlet surface mesh file name
  std::vector<std::string> sur_file_out, sur_file_out_write;
  sur_file_out.resize( num_outlet );
  sur_file_out_write.resize( num_outlet );
  
  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::ostringstream ss, sw;
    ss<<sur_file_out_base;
    sw<<out_out_base;
    if( ii/10 == 0 ) 
    {
      ss<<"00";
      sw<<"00";
    }
    else if( ii/100 == 0 )
    {
      ss<<"0";
      sw<<"0";
    }
    ss<<ii<<".vtp";
    sw<<ii<<".vtp";
    sur_file_out[ii] = ss.str();
    sur_file_out_write[ii] = sw.str();
  }

  // Check the files are on the disk  
  SYS_T::file_check(geo_file);
  
  /*
  SYS_T::file_check(sur_file_in);
  SYS_T::file_check(sur_file_wall);
 
  for( unsigned int ii=0; ii<sur_file_out.size(); ++ii )
    SYS_T::file_check(sur_file_out[ii]);
  */

  std::cout<<"Status: All vtu/vtp files are found.\n";

  // Update the volumetric mesh file
  int nstart, estart;
  SV_T::update_sv_vtu( geo_file, geo_out_name, nstart, estart );

  SYS_T::print_fatal_if( nstart != 1 || estart != 1, "Error: SV file node or element starting index is not 1. Check the files. \n");
  
  std::cout<<"Status: "<<geo_file<< " updated, the starting node index "<<nstart<<" and the starting element index "<<estart<<" are corrected to 0.\n";

  /*
  // Now use the nstart and estart to correct the vtp files  
  SV_T::update_sv_vtp( sur_file_in, inl_out_name, nstart, estart );
  std::cout<<"Status: inflow wall mesh is updated. \n";

  SV_T::update_sv_vtp( sur_file_wall, wal_out_name, nstart, estart );
  std::cout<<"Stauts: wall mesh is updated. \n";

  for( int ii=0; ii<num_outlet; ++ii )
    SV_T::update_sv_vtp( sur_file_out[ii], sur_file_out_write[ii], nstart, estart );
  
  std::cout<<"Status: outflow meshes are updated.\n";

  std::cout<<"Conversion is complete.\n";
  */

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
