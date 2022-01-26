// ============================================================================
// vis_fsi_wss.cpp
//
// This is the visualization driver for WSS.
//
// Date: Jan 26 2022
// ============================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"

int main( int argc, char * argv[] )
{
  int time_index = 0;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");
  
  // Read in the mesh file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string wall_file = cmd_h5r -> read_string("/", "sur_f_file_wall");

  delete cmd_h5r; H5Fclose(prepcmd_file);
  
  hid_t anacmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT );
  
  HDF5_Reader * ana_h5r = new HDF5_Reader( anacmd_file );
  
  const std::string sol_bname = ana_h5r -> read_string("/", "sol_bName");
  
  const double fluid_mu = ana_h5r -> read_doubleScalar("/", "fl_mu");

  delete ana_h5r; H5Fclose(anacmd_file);

  SYS_T::GetOptionInt("-time_index", time_index);

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -time_index: "<<time_index<<endl;
  cout<<"----------------------------------\n";
  cout<<" sol_bname: "<<sol_bname<<endl;
  cout<<" fl_mu: "    <<fluid_mu<<endl;
  cout<<" geo_file: " <<geo_file<<endl;
  cout<<" wall_file: "<<wall_file<<endl;
  cout<<" out_bname: "<<out_bname<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  SYS_T::file_check( geo_file.c_str() );
  SYS_T::file_check( wall_file.c_str() );

  PetscFinalize();
  return EXIT_SUCCESS;
}




// EOF
