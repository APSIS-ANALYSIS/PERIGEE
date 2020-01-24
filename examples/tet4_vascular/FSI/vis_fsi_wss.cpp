// ==================================================================
// vis_fsi_wss.cpp
//
// This is the visualization driver for WSS.
//
// Author: Ju Liu
// Date: Jan 24 2020
// ==================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"


int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  double fluid_mu = 4.0e-2;

  const int dof = 7;

  int time_index = 0;

  std::string geo_file, wall_file;

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1, "ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionInt("-time_index", time_index);

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");

  // Read in the mesh file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  cmd_h5r -> read_string("/", "geo_file", geo_file);
  cmd_h5r -> read_string("/", "sur_f_file_wall", wall_file);

  delete cmd_h5r; H5Fclose(prepcmd_file);

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sol_bname: "<<sol_bname<<endl;
  cout<<" -fl_mu: "<<fluid_mu<<endl;
  cout<<" -time_index: "<<time_index<<endl;
  cout<<"----------------------------------\n";
  cout<<" geo_file: "<<geo_file<<endl;
  cout<<" wall_file: "<<wall_file<<endl;
  cout<<" out_bname: "<<out_bname<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  SYS_T::file_exist_check( geo_file.c_str() );
  SYS_T::file_exist_check( wall_file.c_str() );

  // Read in the whole FSI volumetric mesh  
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN, phy_tag;
  std::vector<double> v_ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), v_nFunc, v_nElem, v_ctrlPts, v_vecIEN, phy_tag);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Readd the wall surface mesh
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN, global_node_idx, global_ele_idx;

  TET_T::read_vtp_grid( wall_file.c_str(), nFunc, nElem, ctrlPts, vecIEN,
      global_node_idx, global_ele_idx );

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";

  // Each surface triangle element requires an additional node: interior_node
  std::vector<int> interior_node; interior_node.resize( nElem );
 
  // Interior nodes's xyz coordinates 
  std::vector<double> interior_node_coord;
  interior_node_coord.resize(3*nElem);

  // Unit outward normal vector
  std::vector< std::vector<double> > outnormal;
  outnormal.resize( nElem );

  // Triangle element surface area
  std::vector<double> tri_area; tri_area.resize( nElem );

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
