// ==================================================================
// This is a file converter driver that reads in the svpre file format
// to locate the volumetric vtu, the inlet vtp, the wall vtp, and
// the outlets' vtp files, and then make their indices strarting from
// zero.
// 
// This drive will read in the files listed in svpre_file.
// The svpre_file will have 1 + 1 + 1 + num_outlet lines
//
// The following is very IMPORTANT!! 
// Make sure the first line is the vtu volume file
// the second line is the exterior wall vtp file
// the third line is the inlet face vtp file
// and the rest lines are the outlet faces vtp files
//
// Author: Ju Liu
// Date: Sept 11 2019
// ==================================================================
#include "SV_Tools.hpp"

int main( int argc, char * argv[] )
{
  std::string svpre_file("");
  std::string geo_file("");
  std::string sur_file_wall("");
  std::string sur_file_in("");
  std::vector<std::string> sur_file_out;

  std::string geo_out_name("whole_vol.vtu");
  std::string wal_out_name("wall_vol.vtp");
  std::string inl_out_name("inflow_vol.vtp");
  std::string out_out_base("outflow_vol_");
  
  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1,"ERROR: converter is a serial routine! \n");
  
  SYS_T::GetOptionString("-svpre_file", svpre_file); 
  SYS_T::GetOptionString("-vol_geo_name", geo_out_name);
  SYS_T::GetOptionString("-inl_geo_name", inl_out_name);
  SYS_T::GetOptionString("-wal_geo_name", wal_out_name);
  SYS_T::GetOptionString("-out_geo_base", out_out_base);

  std::cout<<"==== Command Line Arguments ===="<<std::endl;
  std::cout<<" -svpre_file:   "<<svpre_file<<std::endl;
  std::cout<<"names to be written: \n";
  std::cout<<" -vol_geo_name: "<<geo_out_name<<std::endl;
  std::cout<<" -wal_geo_name: "<<wal_out_name<<std::endl;
  std::cout<<" -inl_geo_name: "<<inl_out_name<<std::endl;
  std::cout<<" -out_geo_base: "<<out_out_base<<std::endl;
  std::cout<<"================================"<<std::endl;

  // Now read in the svpre file
  SYS_T::file_check( svpre_file );

  std::ifstream reader;
  reader.open( svpre_file.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  std::string column_1, column_2;
  
  // Read vtu file
  std::getline(reader, sline);
  sstrm.str( sline );
  sstrm >> column_1;
  sstrm >> column_2;
  sstrm.clear();

  if( column_1.compare("mesh_and_adjncy_vtu") == 0 )
    geo_file = column_2;
  else
    SYS_T::print_fatal("Error: svpre 1st line should be vtu file.\n");

  SYS_T::file_check(geo_file);

  // Read wall vtp file
  std::getline(reader, sline);
  sstrm.str( sline );
  sstrm >> column_1;
  sstrm >> column_2;
  sstrm.clear();

  if( column_1.compare("set_surface_id_vtp") == 0 )
    sur_file_wall = column_2;
  else
    SYS_T::print_fatal("Error: svpre 2nd line should be vtp file.\n");

  SYS_T::file_check(sur_file_wall);
 
  // Read inlet vtp file
  std::getline(reader, sline);
  sstrm.str( sline );
  sstrm >> column_1;
  sstrm >> column_2;
  sstrm.clear();

  if( column_1.compare("set_surface_id_vtp") == 0 )
    sur_file_in = column_2;
  else
    SYS_T::print_fatal("Error: svpre 3rd line should be vtp file.\n");

  SYS_T::file_check(sur_file_in); 
  
  // Read outlet files
  int num_outlet = 0;
  while( std::getline(reader, sline) )
  {
    if( sline.compare("") !=0 )
    {
      sstrm.str( sline );
      sstrm >> column_1;
      sstrm >> column_2;
      sstrm.clear();

      if( column_1.compare("set_surface_id_vtp") == 0 )
      {
        sur_file_out.push_back( column_2 );
        num_outlet += 1;
        SYS_T::file_check(column_2); 
      }
      else
        SYS_T::print_fatal("Error: svpre line should be vtp file.\n");
    }
  }

  reader.close();

  std::cout<<"Status: Finish reading "<<svpre_file<<". All vtu/vtp files are found, and the number of outlet faces is "<<num_outlet<<std::endl;

  // Generate the outlet surface mesh output file name
  std::vector<std::string> sur_file_out_write;
  sur_file_out_write.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::ostringstream sw;
    sw<<out_out_base;
    if( ii/10 == 0 ) sw<<"00";
    else if( ii/100 == 0 ) sw<<"0";

    sw<<ii<<".vtp";
    sur_file_out_write[ii] = sw.str();
  }

  // Update teh volumetric mesh file
  int nstart, estart;
  SV_T::update_sv_vtu( geo_file, geo_out_name, nstart, estart );

  SYS_T::print_fatal_if( nstart != 1 || estart != 1, "Error: SV file node or element starting index is not 1. Check the files.\n");

  std::cout<<"Status: "<<geo_file<< " updated, the starting node index "<<nstart<<" and the starting element index "<<estart<<" are corrected to 0.\n";

  // Now use the nstart and estart to correct the vtp files
  SV_T::update_sv_vtp( sur_file_in, inl_out_name, nstart, estart );
  std::cout<<"Status: inflow wall mesh is updated.\n";

  SV_T::update_sv_vtp( sur_file_wall, wal_out_name, nstart, estart );
  std::cout<<"Stauts: wall mesh is updated.\n";

  for(int ii=0; ii<num_outlet; ++ii)
    SV_T::update_sv_vtp( sur_file_out[ii], sur_file_out_write[ii], nstart, estart );

  std::cout<<"Status: "<<num_outlet<<" outflow face meshes are updated.\n";

  std::cout<<"Conversion is complete.\n";

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
