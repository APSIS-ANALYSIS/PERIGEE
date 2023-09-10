// ==================================================================
// This is a file convert function that correct FSI mesh files node and 
// element indices to make them start from zero. The default starting
// index in SimVascular is 1. This will also adjust the structural
// mesh indices based on the fluid mesh. We default that the fluid
// mesh start with 0 and the solid mesh indices follow those of the 
// fluid mesh.
//
// Author: Ju Liu
// Date: April 03 2019
// ==================================================================
#include "SV_Tools.hpp"

int main( int argc, char * argv[] )
{
  // FSI whole mesh file name
  std::string geo_out_name("whole_vol.vtu");

  // We assume the files contain a vtu file for the volume mesh
  // a single inlet and a single exterior wall mesh in vtp format
  // and several outlet mesh in vtp format, named as: outlet_xxx.vtp
  // Users are allowed to input the geo_file, sur_file_in, and
  // sur_file_wall, as well as sur_file_out_base. 
  std::string geo_f_file("lumen.vtu");
  std::string sur_f_file_wall("lumen_wall.vtp");
  std::string sur_f_file_in_base( "lumen_inlet_");
  std::string sur_f_file_out_base("lumen_outlet_");
  int num_inlet = 1;
  int num_outlet = 1;

  std::string geo_f_out_name("lumen_vol.vtu");
  std::string wal_f_out_name("lumen_wall_vol.vtp");
  std::string inl_f_out_base("lumen_inlet_vol_");
  std::string out_f_out_base("lumen_outlet_vol_");

  std::string geo_s_file("tissue.vtu");
  std::string sur_s_file_wall("tissue_wall.vtp");
  std::string sur_s_file_inner_wall("tissue_inner.vtp");
  std::string sur_s_file_in_base("tissue_inlet_");
  std::string sur_s_file_out_base("tissue_outlet_");

  std::string geo_s_out_name("tissue_vol.vtu");
  std::string wal_s_out_name("tissue_wall_vol.vtp");
  std::string inn_s_out_name("tissue_interior_wall_vol.vtp");
  std::string inl_s_out_base("tissue_inlet_vol_");
  std::string out_s_out_base("tissue_outlet_vol_");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  // Make sure this is a serial run
  SYS_T::print_fatal_if( SYS_T::get_MPI_size() != 1, "ERROR: converter is a serial routine! \n" );
  
  SYS_T::GetOptionInt(   "-num_outlet",     num_outlet);
  SYS_T::GetOptionInt(   "-num_inlet",      num_inlet);
  
  SYS_T::GetOptionString("-vol_geo_name", geo_out_name);
  
  SYS_T::GetOptionString("-vol_f_mesh",     geo_f_file);
  SYS_T::GetOptionString("-wall_f_mesh",    sur_f_file_wall);
  SYS_T::GetOptionString("-inlet_f_mesh",   sur_f_file_in_base);
  SYS_T::GetOptionString("-outlet_f_mesh",  sur_f_file_out_base);
  
  SYS_T::GetOptionString("-vol_s_mesh",     geo_s_file);
  SYS_T::GetOptionString("-wall_s_mesh",    sur_s_file_wall);
  SYS_T::GetOptionString("-inner_s_mesh",   sur_s_file_inner_wall);
  SYS_T::GetOptionString("-inlet_s_mesh",   sur_s_file_in_base);
  SYS_T::GetOptionString("-outlet_s_mesh",  sur_s_file_out_base);
  
  SYS_T::GetOptionString("-vol_geo_f_name", geo_f_out_name);
  SYS_T::GetOptionString("-wal_geo_f_name", wal_f_out_name);
  SYS_T::GetOptionString("-inl_geo_f_name", inl_f_out_base);
  SYS_T::GetOptionString("-out_geo_f_base", out_f_out_base);
  
  SYS_T::GetOptionString("-vol_geo_s_name", geo_s_out_name);
  SYS_T::GetOptionString("-wal_geo_s_name", wal_s_out_name);
  SYS_T::GetOptionString("-inn_geo_s_name", inn_s_out_name);
  SYS_T::GetOptionString("-inl_geo_s_name", inl_s_out_base);
  SYS_T::GetOptionString("-out_geo_s_base", out_s_out_base);
  
  std::cout<<"==== Command Line Arguments ===="<<std::endl;
  std::cout<<" -num_outlet:     "<<num_outlet            <<std::endl;
  std::cout<<" -num_inlet:      "<<num_inlet             <<std::endl;
  std::cout<<" -vol_f_mesh:     "<<geo_f_file            <<std::endl; 
  std::cout<<" -wall_f_mesh:    "<<sur_f_file_wall       <<std::endl; 
  std::cout<<" -inlet_f_mesh:   "<<sur_f_file_in_base    <<std::endl; 
  std::cout<<" -outlet_f_mesh:  "<<sur_f_file_out_base   <<std::endl; 
  std::cout<<" -vol_s_mesh:     "<<geo_s_file            <<std::endl; 
  std::cout<<" -wall_s_mesh:    "<<sur_s_file_wall       <<std::endl; 
  std::cout<<" -inner_s_mesh:   "<<sur_s_file_inner_wall <<std::endl; 
  std::cout<<" -inlet_s_mesh:   "<<sur_s_file_in_base    <<std::endl; 
  std::cout<<" -outlet_s_mesh:  "<<sur_s_file_out_base   <<std::endl; 
  std::cout<<"----- names to be written: \n";
  std::cout<<" -vol_geo_name:   "<<geo_out_name          <<std::endl; 
  std::cout<<" -vol_geo_f_name: "<<geo_f_out_name        <<std::endl;
  std::cout<<" -wal_geo_f_name: "<<wal_f_out_name        <<std::endl;
  std::cout<<" -inl_geo_f_base: "<<inl_f_out_base        <<std::endl;
  std::cout<<" -out_geo_f_base: "<<out_f_out_base        <<std::endl;
  std::cout<<" -vol_geo_s_name: "<<geo_s_out_name        <<std::endl;
  std::cout<<" -wal_geo_s_name: "<<wal_s_out_name        <<std::endl;
  std::cout<<" -inn_geo_s_name: "<<inn_s_out_name        <<std::endl;
  std::cout<<" -inl_geo_s_name: "<<inl_s_out_base        <<std::endl;
  std::cout<<" -out_geo_s_base: "<<out_s_out_base        <<std::endl;
  std::cout<<"================================="         <<std::endl;

  // Boundary check for the number of inlets and outlets
  SYS_T::print_fatal_if(num_inlet  < 1,   "Error: -num_inlet cannot be less than 1.\n");
  SYS_T::print_fatal_if(num_inlet  > 100, "Error: -num_inlet cannot be more than 100.\n");
  
  SYS_T::print_fatal_if(num_outlet < 1,   "Error: -num_outlet cannot be less than 1.\n");
  SYS_T::print_fatal_if(num_outlet > 100, "Error: -num_outlet cannot be more than 100.\n");

  // If the corresponding name are the same, old file will be overwirtten
  // make a check.
  if( geo_f_file == geo_f_out_name || sur_f_file_in_base == inl_f_out_base || sur_f_file_wall == wal_f_out_name || sur_f_file_out_base == out_f_out_base || geo_s_file == geo_s_out_name || sur_s_file_in_base == inl_s_out_base || sur_s_file_wall == wal_s_out_name || sur_s_file_out_base == out_s_out_base ) std::cout<<"Warning: Some files have the same read and write names, new file will replace the old files.\n";

  // Generate the outlet surface mesh file name
  std::vector<std::string> sur_f_file_in( num_inlet ), sur_f_file_in_write( num_inlet );
  std::vector<std::string> sur_s_file_in( num_inlet ), sur_s_file_in_write( num_inlet );
  
  for(int ii=0; ii<num_inlet; ++ii)
  {
    sur_f_file_in[ii]       = SYS_T::gen_capfile_name( sur_f_file_in_base, ii, ".vtp" ); 
    sur_f_file_in_write[ii] = SYS_T::gen_capfile_name( inl_f_out_base,     ii, ".vtp" ); 
    sur_s_file_in[ii]       = SYS_T::gen_capfile_name( sur_s_file_in_base, ii, ".vtp" ); 
    sur_s_file_in_write[ii] = SYS_T::gen_capfile_name( inl_s_out_base,     ii, ".vtp" ); 
  }
  
  std::vector<std::string> sur_f_file_out( num_outlet ), sur_f_file_out_write( num_outlet );
  std::vector<std::string> sur_s_file_out( num_outlet ), sur_s_file_out_write( num_outlet );
  
  for(int ii=0; ii<num_outlet; ++ii)
  {
    sur_f_file_out[ii]       = SYS_T::gen_capfile_name( sur_f_file_out_base, ii, ".vtp" ); 
    sur_f_file_out_write[ii] = SYS_T::gen_capfile_name( out_f_out_base,      ii, ".vtp" ); 
    sur_s_file_out[ii]       = SYS_T::gen_capfile_name( sur_s_file_out_base, ii, ".vtp" ); 
    sur_s_file_out_write[ii] = SYS_T::gen_capfile_name( out_s_out_base,      ii, ".vtp" ); 
  }

  // Check the files are on the disk  
  SYS_T::file_check(geo_f_file);
  SYS_T::file_check(sur_f_file_wall);
 
  SYS_T::file_check(geo_s_file);
  SYS_T::file_check(sur_s_file_wall);
 
  for( int ii=0; ii<num_inlet; ++ii )
  {
    SYS_T::file_check( sur_f_file_in[ii] );
    SYS_T::file_check( sur_s_file_in[ii] );
  }
 
  for( int ii=0; ii<num_outlet; ++ii )
  {
    SYS_T::file_check( sur_f_file_out[ii] );
    SYS_T::file_check( sur_s_file_out[ii] );
  }
  std::cout<<"Status: All vtu/vtp files are found.\n";

  // If there are additional cap file throw an warning
  if( SYS_T::file_exist(SYS_T::gen_capfile_name(sur_f_file_in_base, num_inlet, ".vtp")) ||
      SYS_T::file_exist(SYS_T::gen_capfile_name(sur_s_file_in_base, num_inlet, ".vtu")) )
    cout<<endl<<"Warning: there are additional inlet surface files on disk. Check num_inlet please.\n\n";

  if( SYS_T::file_exist(SYS_T::gen_capfile_name(sur_f_file_out_base, num_outlet, ".vtp")) ||
      SYS_T::file_exist(SYS_T::gen_capfile_name(sur_s_file_out_base, num_outlet, ".vtu")) )
    cout<<endl<<"Warning: there are additional outlet surface files on disk. Check num_outlet please.\n\n";

  const std::vector<int> fluid_global_nodal_id = VTK_T::read_int_PointData( geo_f_file, "GlobalNodeID" );
  for(int ii=0; ii<VEC_T::get_size( fluid_global_nodal_id ); ++ii) SYS_T::print_fatal_if( ii+fluid_global_nodal_id[0] != fluid_global_nodal_id[ii], "Error: fluid global id %d does not match with its natrual id %d. \n", fluid_global_nodal_id[ii], ii );
  
  SYS_T::commPrint("Status: fluid mesh nodal id starts from %d and are listed consecutively.\n", fluid_global_nodal_id[0]);

  const std::vector<int> fluid_global_elem_id = VTK_T::read_int_CellData( geo_f_file, "GlobalElementID" );

  for(int ii=0; ii<VEC_T::get_size( fluid_global_elem_id ); ++ii) SYS_T::print_fatal_if( ii+fluid_global_elem_id[0] != fluid_global_elem_id[ii], "Error: fluid global elem id %d does not match with its natrual id %d. \n", fluid_global_elem_id[ii], ii );

  SYS_T::commPrint("Status: fluid mesh elem id starts from %d and are listed consecutively.\n", fluid_global_elem_id[0]);
  
  const std::vector<int> solid_global_nodal_id = VTK_T::read_int_PointData( geo_s_file, "GlobalNodeID" );
  for(int ii=0; ii<VEC_T::get_size( solid_global_nodal_id ); ++ii) SYS_T::print_fatal_if( ii+solid_global_nodal_id[0] != solid_global_nodal_id[ii], "Error: solid global id %d does not match with its natrual id %d. \n", solid_global_nodal_id[ii], ii);

  SYS_T::commPrint("Status: solid mesh nodal id starts from %d and are listed consecutively.\n", solid_global_nodal_id[0]);
  
  const std::vector<int> solid_global_elem_id = VTK_T::read_int_CellData( geo_s_file, "GlobalElementID" );

  for(int ii=0; ii<VEC_T::get_size( solid_global_elem_id ); ++ii) SYS_T::print_fatal_if( ii+solid_global_elem_id[0] != solid_global_elem_id[ii], "Error: solid global elem id %d does not match with its natrual id %d. \n", solid_global_elem_id[ii], ii );

  SYS_T::commPrint("Status: solid mesh elem id starts from %d and are listed consecutively.\n", solid_global_elem_id[0]);
  
  // Obtain the nstart and estart, the first indices for element and nodal
  // indices. Make sure they start from 1.
  // Write the updated fluid domain mesh into a vtu file.
  int nstart, estart;
  SV_T::update_sv_vtu( geo_f_file, geo_f_out_name, nstart, estart );

  if( nstart == 0 ) std::cout<<"Status: "<<geo_f_file<<" nodal index starts from 0. \n";
  else if( nstart == 1 ) std::cout<<"Status: "<<geo_f_file<<" nodal index starts from 1 and has been updated. \n";
  else std::cout<<"Warning: "<<geo_f_file<<" nodal index starts from "<<nstart<<" and has bee updated. \n";

  if( estart == 0 ) std::cout<<"Status: "<<geo_f_file<<" element index starts from 0. \n";
  else if( estart == 1 ) std::cout<<"Status: "<<geo_f_file<<" element index starts from 1 and has been updated. \n";
  else std::cout<<"Warning: "<<geo_f_file<<" element index starts from "<<estart<<" and has bee updated. \n";

  std::cout<<"Status: "<<geo_f_file<< " is updated to "<<geo_f_out_name<<'\n';

  // Now use the nstart = 1 and estart = 1 to correct the fluid vtp files  
  SV_T::update_sv_vtp( sur_f_file_wall, wal_f_out_name, nstart, estart );
  std::cout<<"Stauts: "<<sur_f_file_wall<<" is updated to "<<wal_f_out_name<<'\n';

  for( int ii=0; ii<num_inlet; ++ii )
  {
    SV_T::update_sv_vtp( sur_f_file_in[ii], sur_f_file_in_write[ii], nstart, estart );
    std::cout<<"Status: "<<sur_f_file_in[ii]<<" is updated to "<<sur_f_file_in_write[ii]<<'\n';
  }

  for( int ii=0; ii<num_outlet; ++ii )
  {
    SV_T::update_sv_vtp( sur_f_file_out[ii], sur_f_file_out_write[ii], nstart, estart );
    std::cout<<"Status: "<<sur_f_file_out[ii]<<" is updated to "<<sur_f_file_out_write[ii]<<'\n';
  }

  // Check the interface meshes are compatible.
  SV_T::compare_sv_vtp(sur_f_file_wall, sur_s_file_inner_wall);

  // Merge the fluid and solid domain into one FSI domain
  std::vector<int> map_s_node, map_s_elem;
  SV_T::gen_sv_fsi_vtus( geo_f_file, geo_s_file, wal_f_out_name, geo_out_name,
      geo_s_out_name, map_s_node, map_s_elem );

  // Update the solid vtp files using the mappings generated.
  SV_T::update_sv_vtp( sur_s_file_wall, wal_s_out_name, nstart, estart, map_s_node, map_s_elem );
  std::cout<<"Stauts: "<<sur_s_file_wall<<" is updated to "<<wal_s_out_name<<'\n';

  SV_T::update_sv_vtp( sur_s_file_inner_wall, inn_s_out_name, nstart, estart, map_s_node, map_s_elem );
  std::cout<<"Stauts: "<<sur_s_file_inner_wall<<" is updated to "<<inn_s_out_name<<'\n';

  for( int ii=0; ii<num_inlet; ++ii )
  {
    SV_T::update_sv_vtp( sur_s_file_in[ii], sur_s_file_in_write[ii], nstart, estart, map_s_node, map_s_elem );
    std::cout<<"Status: "<<sur_s_file_in[ii]<<" is updated to "<<sur_s_file_in_write[ii]<<'\n';
  }

  for( int ii=0; ii<num_outlet; ++ii )
  {
    SV_T::update_sv_vtp( sur_s_file_out[ii], sur_s_file_out_write[ii], nstart, estart, map_s_node, map_s_elem );
    std::cout<<"Status: "<<sur_s_file_out[ii]<<" is updated to "<<sur_s_file_out_write[ii]<<'\n';
  }
  
  // Finalize this routine
  std::cout<<"Conversion is complete. \n";
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
