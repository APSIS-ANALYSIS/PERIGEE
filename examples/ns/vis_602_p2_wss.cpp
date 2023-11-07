// ==================================================================
// vis_602_p2_wss.cpp
//
// WSS visualization for 27-node hex elements.
//
// Date: Oct 30 2023
// ==================================================================
#include "Hex_Tools.hpp"
#include "QuadPts_vis_quad9.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_vis_hex27.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"

std::vector<int> range_generator( const int &ii, const int &jj, const int &kk, const int &ll );

std::vector<int> ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size );

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof );

int get_quad_local_id( const double * const &coor_x,
    const double * const &coor_y,
    const double * const &coor_z,
    const int &len,
    const double &x, const double &y, const double &z,
    const double &tol = 1.0e-8 );

void write_quad_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< Vector_3 > &wss_on_node );

void write_quad_grid_tawss_osi( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<double> &tawss,
    const std::vector<double> &osi );

int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  constexpr int nLocBas = 9;
  constexpr int v_nLocBas = 27;

  constexpr int dof = 4; 

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::print_fatal_if( SYS_T::get_MPI_size() != 1, "ERROR: vis_602_p2_wss is a serial program! \n");

  // Read the geometry file name from preprocessor hdf5 file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  const std::string geo_file  = cmd_h5r -> read_string("/", "geo_file");
  const std::string wall_file = cmd_h5r -> read_string("/", "sur_file_wall");
  const int elemType = cmd_h5r -> read_intScalar("/", "elemType");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Read the material property from the solver HDF5 file
  prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  cmd_h5r = new HDF5_Reader( prepcmd_file );

  const double fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Enforce the element to be triquadratic hex for now
  if( elemType != 602 ) SYS_T::print_fatal("Error: element type should be 602 triquadratic hex element.\n");

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
  SYS_T::file_check( geo_file );
  SYS_T::file_check( wall_file );

  // Read the volumetric mesh
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Read the wall surface mesh
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN;

  VTK_T::read_vtu_grid(wall_file, nFunc, nElem, ctrlPts, vecIEN);
  
  const std::vector<int> global_node_idx = VTK_T::read_int_PointData( wall_file, "GlobalNodeID");
  const std::vector<int> global_ele_idx = VTK_T::read_int_CellData( wall_file, "GlobalElementID");

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";

  std::vector<int> interior_node{}; 

  std::vector<double> interior_node_coord(3*4*nElem, 0.0); // 4 means that each wall surface element has 4 interior nodes

  std::vector<int> interior_node_local_index{};

  for(int ee=0; ee<nElem; ++ee)
  {
    std::vector<int> quadn { global_node_idx[ vecIEN[nLocBas*ee+0] ], global_node_idx[ vecIEN[nLocBas*ee+1] ], 
                             global_node_idx[ vecIEN[nLocBas*ee+2] ], global_node_idx[ vecIEN[nLocBas*ee+3] ] };
    
    const int hexn[8] { v_vecIEN[ global_ele_idx[ee]*v_nLocBas+0 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+1 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+2 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+3 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+4 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+5 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+6 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+7 ] };

    int node_check = 0;
    for(int ii=0; ii<8; ++ii)
    {
      const bool gotnode = VEC_T::is_invec( quadn, hexn[ii] );
  
      if(!gotnode)
      { 
        interior_node.push_back(hexn[ii]); // interior node's global volumetric mesh nodal index
        interior_node_local_index.push_back(ii); // the local indices of nodes on the wall surface
      }
      else 
        node_check += 1;
    }

    SYS_T::print_fatal_if(node_check != 4, "Error: the associated hex element is incompatible with the quad element.\n");
    
    // Now we have found the interior node's volumetric mesh index, record its
    // spatial xyz coordinate
    for (int ii=0; ii<4; ++ii)
    {
      for (int jj=0; jj<3; ++jj)
      {
        interior_node_coord[3*(4*ee + ii )+jj] = v_ctrlPts[ 3*interior_node[4*ee + ii] + jj ];
      }
    } 
  }

  SYS_T::print_fatal_if(VEC_T::get_size(interior_node) != 4*nElem, "Error: the length of the interior_node vector is incorrect.\n");
  SYS_T::print_fatal_if(VEC_T::get_size(interior_node_local_index) != 4*nElem, "Error: the length of the interior_node_local_index vector is incorrect.\n");

  // Volumetric element visualization sampling point 
  IQuadPts * quad = new QuadPts_vis_hex27();

  quad -> print_info();

  FEAElement * element = new FEAElement_Hex27( quad-> get_num_quadPts() );

  double * v_ectrl_x = new double [v_nLocBas];
  double * v_ectrl_y = new double [v_nLocBas];
  double * v_ectrl_z = new double [v_nLocBas];
  double * esol_u  = new double [v_nLocBas];
  double * esol_v  = new double [v_nLocBas];
  double * esol_w  = new double [v_nLocBas];

  IQuadPts * quad_vis = new QuadPts_vis_quad9();

  quad_vis -> print_info();

  IQuadPts * quad_gau = new QuadPts_Gauss_Quad( quad_vis->get_num_quadPts_x(), quad_vis->get_num_quadPts_y() );

  quad_gau -> print_info();

  FEAElement * element_quad = new FEAElement_Quad9_3D_der0( quad_vis->get_num_quadPts() );

  // Read the mappings of the nodal indices
  const std::vector<int> analysis_new2old = ReadNodeMapping("node_mapping.h5", "new_2_old", v_nFunc );

  double * Rx = new double [v_nLocBas];
  double * Ry = new double [v_nLocBas];
  double * Rz = new double [v_nLocBas];

  // Container for Time averaged WSS and OSI
  std::vector<double> tawss( nFunc, 0.0 ); 
  std::vector<double> osi( nFunc, 0.0 ); 
  std::vector< Vector_3 > osi_top( nFunc, Vector_3(0.0, 0.0, 0.0) );

  const double inv_T = 1.0 / ( static_cast<double>((time_end - time_start)/time_step) + 1.0 );

  // Loop over time
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    std::ostringstream time_index;
    time_index.str("");
    time_index << 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s and Write %s \n", time, name_to_read.c_str(), name_to_write.c_str());

    // Read the solution vector and renumber them based on the nodal mappings
    const auto sol = ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof );

    // Container for (averaged) WSS
    std::vector< Vector_3 > wss_ave( nFunc, Vector_3(0.0, 0.0, 0.0) );

    // Container for the element area associated with surface nodes
    std::vector<double> node_area( nFunc, 0.0 );

    for(int ee=0; ee<nElem; ++ee)
    {
      // ee element's volumetric element id
      const int ee_vol_id = global_ele_idx[ ee ];

      // get the hexahedron's control points and associated solution nodal
      // value
      for(int ii=0; ii<v_nLocBas; ++ii)
      {
        v_ectrl_x[ii] = v_ctrlPts[ 3 * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 0 ];
        v_ectrl_y[ii] = v_ctrlPts[ 3 * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 1 ];
        v_ectrl_z[ii] = v_ctrlPts[ 3 * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 2 ];

        esol_u[ii] = sol[ dof * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 1 ];
        esol_v[ii] = sol[ dof * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 2 ];
        esol_w[ii] = sol[ dof * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 3 ];
      }

      // Construct the triquadratic hexahedron element
      element -> buildBasis(quad, v_ectrl_x, v_ectrl_y, v_ectrl_z);

      // Obtain the local indices of nodes on the wall surface based on the local indices of the four interior nodes
      const std::vector<int> id_range = range_generator( interior_node_local_index[4*ee], interior_node_local_index[4*ee + 1], interior_node_local_index[4*ee + 2], interior_node_local_index[4*ee + 3] );

      // Obtain the control point coordinates for this element
      double * ectrl_x = new double [nLocBas];
      double * ectrl_y = new double [nLocBas];
      double * ectrl_z = new double [nLocBas];

      // Here the coordinates of the control points on the wall are obtained via the v_vecIEN
      // Note that their order is determined by the id_range
      for(int ii=0; ii<nLocBas; ++ii)
      {
        ectrl_x[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 0 ];
        ectrl_y[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 1 ];
        ectrl_z[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 2 ];
      }

      // Build a basis based on the visualization sampling point for wall
      // quad element
      element_quad -> buildBasis(quad_vis, ectrl_x, ectrl_y, ectrl_z);

      std::vector< Vector_3 > outnormal( nLocBas, Vector_3(0.0, 0.0, 0.0) );

      // For each nodal point, calculate the outward normal vector using the
      // biquadratic quad element. The quad element's ii-th basis corresponds to the
      // hex element's id_range[ii]-th basis
      for(int ii=0; ii<nLocBas; ++ii)
      {
        double len;
        const Vector_3 sur_pt( ectrl_x[ii], ectrl_y[ii], ectrl_z[ii] );
        const Vector_3 int_pt( interior_node_coord[3*(4*ee)+0], interior_node_coord[3*(4*ee)+1], interior_node_coord[3*(4*ee)+2] );

        // id_range[ii] 's outward normal
        outnormal[ii] = element_quad -> get_normal_out( ii, sur_pt, int_pt, len );
      }

      // Now calcualte the element surface area
      element_quad -> buildBasis( quad_gau, ectrl_x, ectrl_y, ectrl_z );

      double quad_area = 0.0;
      for(int qua=0; qua<quad_gau->get_num_quadPts(); ++qua)
        quad_area += element_quad->get_detJac(qua) * quad_gau->get_qw(qua);

      // Here the coordinates of the control points read through the surface-IEN 
      // are just used to find out their global number on the wall
      for(int ii=0; ii<nLocBas; ++ii)
      {
        ectrl_x[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 0 ];
        ectrl_y[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 1 ];
        ectrl_z[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 2 ];
      }

      // Loop over the 9 sampling points on the wall biquadratic quadrangle element
      for(int qua=0; qua<nLocBas; ++qua)
      {
        // Obtain the 27 basis function's value at the wall boundary points
        element -> get_gradR( id_range[qua], Rx, Ry, Rz );

        double ux = 0.0, uy = 0.0, uz = 0.0;
        double vx = 0.0, vy = 0.0, vz = 0.0;
        double wx = 0.0, wy = 0.0, wz = 0.0;

        for(int ii=0; ii<v_nLocBas; ++ii)
        {
          ux += esol_u[ii] * Rx[ii];
          uy += esol_u[ii] * Ry[ii];
          uz += esol_u[ii] * Rz[ii];

          vx += esol_v[ii] * Rx[ii];
          vy += esol_v[ii] * Ry[ii];
          vz += esol_v[ii] * Rz[ii];

          wx += esol_w[ii] * Rx[ii];
          wy += esol_w[ii] * Ry[ii];
          wz += esol_w[ii] * Rz[ii];
        }

        // obtain the hex element's id_range[qua] node's outward normal vector 
        const double nx = outnormal[qua].x();
        const double ny = outnormal[qua].y();
        const double nz = outnormal[qua].z();

        const double ax = 2.0 * ux * nx + (uy + vx) * ny + (uz + wx) * nz;
        const double ay = (vx + uy) * nx + 2.0 * vy * ny + (vz + wy) * nz;
        const double az = (wx + uz) * nx + (wy + vz) * ny + 2.0 * wz * nz;

        const double b = ax * nx + ay * ny + az * nz;

        const double wss_x = fluid_mu * ( ax - b * nx );
        const double wss_y = fluid_mu * ( ay - b * ny );
        const double wss_z = fluid_mu * ( az - b * nz );

        const int quad_local_id = get_quad_local_id( ectrl_x, ectrl_y, ectrl_z,
            nLocBas, v_ectrl_x[ id_range[qua] ], v_ectrl_y[ id_range[qua] ], 
            v_ectrl_z[ id_range[qua] ], 1.0e-8 );
       
        const int quad_global_id = vecIEN[nLocBas * ee + quad_local_id];

        wss_ave[ quad_global_id ]   += quad_area * Vector_3(wss_x, wss_y, wss_z);
        node_area[ quad_global_id ] += quad_area;

      } // loop over the sampling points (on surface)

      delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;

    } // End of loop over element

    for(int ii=0; ii<nFunc; ++ii) wss_ave[ii] *= (1.0 / node_area[ii]);

    // Write the wall shear stress at this time instance
    write_quad_grid_wss( name_to_write, nFunc, nElem, ctrlPts, vecIEN, wss_ave );

    for(int ii=0; ii<nFunc; ++ii)
    {
      tawss[ii]   += inv_T * wss_ave[ii].norm2();
      osi_top[ii] += inv_T * wss_ave[ii];
    }

  } // End of loop over time

  MPI_Barrier(PETSC_COMM_WORLD);

  for(int ii=0; ii<nFunc; ++ii)
  {
    const double mag = osi_top[ii].norm2();

    // We will disallow very small tawss in the osi calculation
    if( std::abs(tawss[ii]) > 1.0e-12 )
      osi[ii] = 0.5 * (1.0 - mag / tawss[ii] );
    else
      osi[ii] = 0.0;
  }

  // write the TAWSS and OSI
  std::string tawss_osi_file("SOL_TAWSS_OSI" );
  write_quad_grid_tawss_osi( tawss_osi_file, nFunc, nElem, ctrlPts, vecIEN, tawss, osi );

  delete [] v_ectrl_x; delete [] v_ectrl_y; delete [] v_ectrl_z;
  delete [] esol_u; delete [] esol_v; delete [] esol_w;
  delete [] Rx; delete [] Ry; delete [] Rz; 
  delete quad; delete element;
  delete quad_vis; delete quad_gau; delete element_quad;

  PetscFinalize();
  return EXIT_SUCCESS;
}

std::vector<int> range_generator( const int &ii, const int &jj, const int &kk, const int &ll )
{
  std::vector<int> surface_id_range(9, -1);

  int interior_id_range[4] {ii, jj, kk, ll};

  std::sort(interior_id_range, interior_id_range + 4);

  const int zeroth[4] {0,1,2,3};
  const int first [4] {4,5,6,7};
  const int second[4] {0,1,4,5};
  const int third [4] {1,2,5,6};
  const int fourth[4] {2,3,6,7};
  const int fifth [4] {0,3,4,7};

  if (std::equal(interior_id_range, interior_id_range + 4, zeroth))
  {
    surface_id_range[0] = 4;
    surface_id_range[1] = 5;
    surface_id_range[2] = 6;
    surface_id_range[3] = 7;
    surface_id_range[4] = 12;
    surface_id_range[5] = 13;
    surface_id_range[6] = 14;
    surface_id_range[7] = 15;
    surface_id_range[8] = 25;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, first))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 3;
    surface_id_range[2] = 2;
    surface_id_range[3] = 1;
    surface_id_range[4] = 11;
    surface_id_range[5] = 10;
    surface_id_range[6] = 9;
    surface_id_range[7] = 8;
    surface_id_range[8] = 24;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, second))
  {
    surface_id_range[0] = 2;
    surface_id_range[1] = 3;
    surface_id_range[2] = 7;
    surface_id_range[3] = 6;
    surface_id_range[4] = 10;
    surface_id_range[5] = 19;
    surface_id_range[6] = 14;
    surface_id_range[7] = 18;
    surface_id_range[8] = 23;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, third))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 4;
    surface_id_range[2] = 7;
    surface_id_range[3] = 3;
    surface_id_range[4] = 16;
    surface_id_range[5] = 15;
    surface_id_range[6] = 19;
    surface_id_range[7] = 11;
    surface_id_range[8] = 20;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, fourth))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 1;
    surface_id_range[2] = 5;
    surface_id_range[3] = 4;
    surface_id_range[4] = 8;
    surface_id_range[5] = 17;
    surface_id_range[6] = 12;
    surface_id_range[7] = 16;
    surface_id_range[8] = 22;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, fifth))
  {
    surface_id_range[0] = 1;
    surface_id_range[1] = 2;
    surface_id_range[2] = 6;
    surface_id_range[3] = 5;
    surface_id_range[4] = 9;
    surface_id_range[5] = 18;
    surface_id_range[6] = 13;
    surface_id_range[7] = 17;
    surface_id_range[8] = 21;
  }
  else
    SYS_T::print_fatal("Error: the interior node index is wrong!\n");

  return surface_id_range;
}

std::vector<int> ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size )
{
  hid_t file_id = H5Fopen(node_mapping_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT);

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    SYS_T::commPrint("Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0];

  if( int(dSize) != node_size )
  {
    SYS_T::commPrint("Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<int> out(node_size, -1);

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, &out[0] );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);

  return out;
}

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof )
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(),
      FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the solution length
  PetscInt get_sol_temp_size;
  VecGetSize(sol_temp, &get_sol_temp_size);
  if( get_sol_temp_size != vec_size )
  {
    SYS_T::commPrint("The solution size %d is not compatible with the size %d given by partition file! \n",
        get_sol_temp_size, vec_size);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<double> veccopy(vec_size, 0.0);
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii)
    veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // copy the solution varibles to the correct location
  std::vector<double> sol(vec_size, 0.0);

  // check the nodemap size
  if( (int)nodemap.size() * in_dof != vec_size ) SYS_T::print_fatal("Error: node map size is incompatible with the solution length. \n");

  for(unsigned int ii=0; ii<nodemap.size(); ++ii)
  {
    const int index = nodemap[ii];
    for(int jj=0; jj<in_dof; ++jj)
      sol[in_dof*index+jj] = veccopy[in_dof*ii+jj];
  }

  return sol;
}

int get_quad_local_id( const double * const &coor_x,
    const double * const &coor_y,
    const double * const &coor_z,
    const int &len,
    const double &x, const double &y, const double &z,
    const double &tol )
{
  for(int ii=0; ii<len; ++ii)
  {
    double dist = 0.0;
    dist  += (x - coor_x[ii]) * (x - coor_x[ii]);
    dist  += (y - coor_y[ii]) * (y - coor_y[ii]);
    dist  += (z - coor_z[ii]) * (z - coor_z[ii]);

    dist = std::sqrt(dist);

    if(dist < tol) return ii;
  }

  SYS_T::print_fatal("Error in get_quad_local_id.\n");

  return -1;
}

void write_quad_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< Vector_3 > &wss_on_node )
{
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // generate the biquadratic quad grid
  HEX_T::gen_quadratic_quad_grid(grid_w, numpts, numcels, pt, ien_array);

  // write wss
  VTK_T::add_Vector3_PointData(grid_w, wss_on_node, "WSS");

  // write vtu
  VTK_T::write_vtkPointSet(filename, grid_w, true);

  grid_w->Delete();
}

void write_quad_grid_tawss_osi( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<double> &tawss,
    const std::vector<double> &osi )
{
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // generate the biquadratic_quad grid
  HEX_T::gen_quadratic_quad_grid(grid_w, numpts, numcels, pt, ien_array);

  // write tawss
  VTK_T::add_double_PointData(grid_w, tawss, "TAWSS");

  // write osi
  VTK_T::add_double_PointData(grid_w, osi, "OSI");

  // write vtu
  VTK_T::write_vtkPointSet(filename, grid_w, true);

  grid_w->Delete();
}

// EOF
