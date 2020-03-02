// ==================================================================
// vis_tet10_wss.cpp
//
// WSS visualization for ten-node tet elements.
// ==================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle6_3d_der0.hpp"

int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  std::string geo_file, wall_file;
  int elemType = 502;

  double fluid_mu = 3.5e-2;

  const int dof = 4; 

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1, "ERROR: vis_tet10_wss is a serial program! \n");
  
  // Read the geometry file name from preprocessor hdf5 file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  cmd_h5r -> read_string("/", "geo_file", geo_file);
  cmd_h5r -> read_string("/", "sur_file_wall", wall_file);
  elemType = cmd_h5r -> read_intScalar("/", "elemType");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Read the material property from the solver HDF5 file
  prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  cmd_h5r = new HDF5_Reader( prepcmd_file );

  fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(prepcmd_file);
  
  // Enforce the element to be quadratic tet for now
  if( elemType != 502 ) SYS_T::print_fatal("Error: element type should be 502 quadratic tet element.\n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");
  
  // Print the command line argument on screen
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sol_bname: "<<sol_bname<<endl;
  cout<<" -time_start: "<<time_start<<endl;
  cout<<" -time_step: "<<time_step<<endl;
  cout<<" -time_end: "<<time_end<<endl;
  cout<<"----------------------------------\n";
  cout<<" geo_file: "<<geo_file<<endl;
  cout<<" wall_file: "<<wall_file<<endl;
  cout<<" elemType: "<<elemType<<endl;
  cout<<" out_bname: "<<out_bname<<endl;
  cout<<" fl_mu: "<<fluid_mu<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Make sure the geometry files exist
  SYS_T::file_exist_check( geo_file.c_str() );
  SYS_T::file_exist_check( wall_file.c_str() );

  // Read the volumetric mesh
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Read the wall surface mesh
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN, global_node_idx, global_ele_idx;

  TET_T::read_vtu_grid( wall_file.c_str(), nFunc, nElem, ctrlPts, vecIEN,
      global_node_idx, global_ele_idx );

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";

  // Now calculate the outward normal vector and element area
  std::vector<int> interior_node;
  interior_node.resize( nElem );

  std::vector<double> interior_node_coord;
  interior_node_coord.resize(3*nElem);

  std::vector< std::vector<double> > outnormal;
  outnormal.resize( 6*nElem ); // quad triangle element has 6 nodes with different normal vectors

  std::vector<double> tri_area;
  tri_area.resize( nElem );

  for(int ee=0; ee<nElem; ++ee)
  {
    std::vector<int> trn; trn.resize(3);
    int ten[4];

    trn[0] = global_node_idx[ vecIEN[3*ee+0] ];
    trn[1] = global_node_idx[ vecIEN[3*ee+1] ];
    trn[2] = global_node_idx[ vecIEN[3*ee+2] ];

    ten[0] = v_vecIEN[ global_ele_idx[ee]*4+0 ];
    ten[1] = v_vecIEN[ global_ele_idx[ee]*4+1 ];
    ten[2] = v_vecIEN[ global_ele_idx[ee]*4+2 ];
    ten[3] = v_vecIEN[ global_ele_idx[ee]*4+3 ];
    
    bool gotnode[4];
    int node_check = 0;
    for(int ii=0; ii<4; ++ii)
    {
      gotnode[ii] = VEC_T::is_invec( trn, ten[ii] );

      if(!gotnode[ii]) interior_node[ee] = ten[ii];
      else node_check += 1;
    }

    SYS_T::print_fatal_if(node_check!=3, "Error: the associated tet element is incompatible with the triangle element.\n");
    
    // Now we have found the interior node's volumetric mesh index, record its
    // coordinate
    interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
    interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
    interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];
 
     
  
  
  
  
  
  
  
  }






























  PetscFinalize();
  return EXIT_SUCCESS;
}




// EOF
