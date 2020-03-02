// ==================================================================
// vis_tet10_wss.cpp
//
// WSS visualization for ten-node tet elements.
// ==================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tri6.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"

void range_generator( const int &ii, std::vector<int> &surface_id_range );


int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  std::string geo_file, wall_file;
  int elemType = 502;

  const int nLocBas = 6;
  const int v_nLocBas = 10;

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
  //if( elemType != 502 ) SYS_T::print_fatal("Error: element type should be 502 quadratic tet element.\n");

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

  // take values 0, 1, 2, or 3, it determines the side of the surface in the
  // tetrahedron
  std::vector<int> interior_node_local_index; 
  interior_node_local_index.resize( nElem );

  std::vector< std::vector<double> > outnormal;
  outnormal.resize( nLocBas*nElem ); // quad triangle element has 6 nodes with different normal vectors

  std::vector<double> tri_area;
  tri_area.resize( nElem );

  IQuadPts * quad_tri_vis = new QuadPts_vis_tri6();

  quad_tri_vis -> print_info();

  IQuadPts * quad_tri_gau = new QuadPts_Gauss_Triangle( quad_tri_vis->get_num_quadPts() );

  quad_tri_gau -> print_info();

  FEAElement * element_tri = new FEAElement_Triangle6_3D_der0( quad_tri_vis-> get_num_quadPts() );

  for(int ee=0; ee<nElem; ++ee)
  {
    std::vector<int> trn; trn.resize(3);
    int ten[4];

    trn[0] = global_node_idx[ vecIEN[nLocBas*ee+0] ];
    trn[1] = global_node_idx[ vecIEN[nLocBas*ee+1] ];
    trn[2] = global_node_idx[ vecIEN[nLocBas*ee+2] ];

    ten[0] = v_vecIEN[ global_ele_idx[ee]*v_nLocBas+0 ];
    ten[1] = v_vecIEN[ global_ele_idx[ee]*v_nLocBas+1 ];
    ten[2] = v_vecIEN[ global_ele_idx[ee]*v_nLocBas+2 ];
    ten[3] = v_vecIEN[ global_ele_idx[ee]*v_nLocBas+3 ];
    
    bool gotnode[4];
    int node_check = 0;
    for(int ii=0; ii<4; ++ii)
    {
      gotnode[ii] = VEC_T::is_invec( trn, ten[ii] );

      if(!gotnode[ii]) 
      {
        interior_node[ee] = ten[ii]; // interior node's global volumetric mesh nodal index
        interior_node_local_index[ee] = ii; // interior node's local tetrahedral element index
      }
      else node_check += 1;
    }

    SYS_T::print_fatal_if(node_check!=3, "Error: the associated tet element is incompatible with the triangle element.\n");
    
    // Now we have found the interior node's volumetric mesh index, record its
    // coordinate
    interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
    interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
    interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];

    // Obtain the control point coordinates for this element
    double * ectrl_x = new double [nLocBas];
    double * ectrl_y = new double [nLocBas];
    double * ectrl_z = new double [nLocBas];

    for(int ii=0; ii<nLocBas; ++ii)
    {
      ectrl_x[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii] + 0 ];
      ectrl_y[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii] + 1 ];
      ectrl_z[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii] + 2 ];
    }

    // Build a basis based on the visualization sampling point 
    element_tri -> buildBasis(quad_tri_vis, ectrl_x, ectrl_y, ectrl_z); 

    // For each nodal point, calculate the outward normal vector 
    for(int ii=0; ii<nLocBas; ++ii)
    {
      double nx, ny, nz, len;

      element_tri -> get_normal_out( ii, ectrl_x[ii], ectrl_y[ii], ectrl_z[ii],
          interior_node_coord[3*ee+0], interior_node_coord[3*ee+1], 
          interior_node_coord[3*ee+2], nx, ny, nz, len );

      outnormal[nLocBas * ee + ii].resize(3);
      outnormal[nLocBas * ee + ii][0] = nx;
      outnormal[nLocBas * ee + ii][1] = ny;
      outnormal[nLocBas * ee + ii][2] = nz;
    } 
  
    // Now calcualte the element surface area
    element_tri -> buildBasis( quad_tri_gau, ectrl_x, ectrl_y, ectrl_z );

    tri_area[ee] = 0.0;
    for(int qua=0; qua<quad_tri_gau->get_num_quadPts(); ++qua)
      tri_area[ee] += element_tri->get_detJac(qua) * quad_tri_gau->get_qw(qua);
 
    delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z; 
  }

  delete quad_tri_vis; delete quad_tri_gau; delete element_tri;
  
  // Volumetric element visualization sampling point 
  IQuadPts * quad = new QuadPts_vis_tet10_v2();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet10_v2( quad-> get_num_quadPts() );

  double * ectrl_x = new double [v_nLocBas];
  double * ectrl_y = new double [v_nLocBas];
  double * ectrl_z = new double [v_nLocBas];
  double * esol_u  = new double [v_nLocBas];
  double * esol_v  = new double [v_nLocBas];
  double * esol_w  = new double [v_nLocBas];



















  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  delete [] esol_u; delete [] esol_v; delete [] esol_w;
  delete quad; delete element;

  PetscFinalize();
  return EXIT_SUCCESS;
}


void range_generator( const int &ii, std::vector<int> &surface_id_range )
{
  surface_id_range.resize(6);
  switch (ii)
  {
    case 0:
      surface_id_range[0] = 1;
      surface_id_range[1] = 2;
      surface_id_range[2] = 3;
      surface_id_range[3] = 5;
      surface_id_range[4] = 8;
      surface_id_range[5] = 9;
      break;
    case 1:
      surface_id_range[0] = 0;
      surface_id_range[1] = 2;
      surface_id_range[2] = 3;
      surface_id_range[3] = 6;
      surface_id_range[4] = 7;
      surface_id_range[5] = 9;
      break;
    case 2:
      surface_id_range[0] = 0;
      surface_id_range[1] = 1;
      surface_id_range[2] = 3;
      surface_id_range[3] = 4;
      surface_id_range[4] = 7;
      surface_id_range[5] = 8;
      break;
    case 3:
      surface_id_range[0] = 0;
      surface_id_range[1] = 1;
      surface_id_range[2] = 2;
      surface_id_range[3] = 4;
      surface_id_range[4] = 5;
      surface_id_range[5] = 6;
      break;
    default:
      SYS_T::print_fatal("Error: the interior node index is wrong!\n");
      break;
  }
}

// EOF
