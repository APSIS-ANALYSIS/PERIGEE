// ==================================================================
// vis_fsi_wss_hex8.cpp
//
// WSS visualization for 8-node hex elements.
//
// Date: Nov 7 2023
// ==================================================================
#include "HDF5_Tools.hpp"
#include "Hex_Tools.hpp"
#include "QuadPts_vis_quad4.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_vis_hex8.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"

std::vector<int> range_generator( const int &ii, const int &jj, const int &kk, const int &ll );

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
  // visualization sampling pattern over time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  constexpr int nLocBas = 4;
  constexpr int v_nLocBas = 8;

  constexpr int dof_v = 3;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");
  
  // Read in the mesh file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string wall_file = cmd_h5r -> read_string("/", "sur_f_file_wall");
  const std::string elemType_str = cmd_h5r -> read_string("/", "elemType");
  const FEType elemType = FE_T::to_FEType(elemType_str);

  delete cmd_h5r; H5Fclose(prepcmd_file);
  
  hid_t anacmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT );
  
  HDF5_Reader * ana_h5r = new HDF5_Reader( anacmd_file );
  
  const std::string sol_bname = ana_h5r -> read_string("/", "sol_bName");
  
  const double fluid_mu = ana_h5r -> read_doubleScalar("/", "fl_mu");

  delete ana_h5r; H5Fclose(anacmd_file);

  // Enforce the element to be trilinear hex for now
  if( elemType != FEType::Hex8 ) SYS_T::print_fatal("Error: element type should be trilinear hex element.\n");

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);  

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -time_start: "<<time_start<<endl;
  cout<<" -time_step: "<<time_step<<endl;
  cout<<" -time_end: "<<time_end<<endl;
  cout<<"----------------------------------\n";
  cout<<" sol_bname: "<<sol_bname<<endl;
  cout<<" fl_mu: "    <<fluid_mu<<endl;
  cout<<" geo_file: " <<geo_file<<endl;
  cout<<" wall_file: "<<wall_file<<endl;
  cout<<" out_bname: "<<out_bname<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  SYS_T::file_check( geo_file );
  SYS_T::file_check( wall_file );

  // ----------------------------------------------------------------
  // Read in the whole FSI volumetric mesh
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);

  const std::vector<int> phy_tag = VTK_T::read_int_CellData(geo_file, "Physics_tag");

  cout<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Read the wall surface mesh
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN;

  VTK_T::read_vtp_grid( wall_file, nFunc, nElem, ctrlPts, vecIEN );

  // They store the coordinates of the control points before deformation
  const std::vector<double> v_ctrlPts_origin(v_ctrlPts);
  const std::vector<double> ctrlPts_origin(ctrlPts);

  const std::vector<int> global_node_idx = VTK_T::read_int_PointData(wall_file, "GlobalNodeID");
  const std::vector<int> global_ele_idx = VTK_T::read_int_CellData(wall_file, "GlobalElementID");

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";

  // Each surface quadrangle element requires 4 additional nodes: interior_node
  std::vector<int> interior_node{};

  // The local indices of nodes on the wall surface  
  std::vector<int> interior_node_local_index{};

  // Identify the interior node for surface elements
  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> quadn { global_node_idx[ vecIEN[nLocBas*ee+0] ], global_node_idx[ vecIEN[nLocBas*ee+1] ],
                                   global_node_idx[ vecIEN[nLocBas*ee+2] ], global_node_idx[ vecIEN[nLocBas*ee+3] ] };

    const int hexn[8] { v_vecIEN[ global_ele_idx[ee]*v_nLocBas+0 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+1 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+2 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+3 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+4 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+5 ],
                        v_vecIEN[ global_ele_idx[ee]*v_nLocBas+6 ], v_vecIEN[ global_ele_idx[ee]*v_nLocBas+7 ] };    

    int node_check = 0;
    for(int ii=0; ii<8; ++ii)
    {
      if( !VEC_T::is_invec( quadn, hexn[ii] ) ) 
      {
      	interior_node.push_back(hexn[ii]);
        interior_node_local_index.push_back(ii);    
      }
      else node_check += 1;
    }

    SYS_T::print_fatal_if(node_check!=4, "Error: the associated hex element is incompatible with the quad element.\n");
  }

  SYS_T::print_fatal_if(VEC_T::get_size(interior_node) != 4*nElem, "Error: the length of the interior_node vector is incorrect.\n");
  SYS_T::print_fatal_if(VEC_T::get_size(interior_node_local_index) != 4*nElem, "Error: the length of the interior_node_local_index vector is incorrect.\n");

  // Volumetric element visualization sampling points and element container
  IQuadPts * quad = new QuadPts_vis_hex8();

  quad -> print_info();

  FEAElement * element = new FEAElement_Hex8( quad-> get_num_quadPts() );

  double * v_ectrl_x = new double [v_nLocBas];
  double * v_ectrl_y = new double [v_nLocBas];
  double * v_ectrl_z = new double [v_nLocBas];
  double * esol_u  = new double [v_nLocBas];
  double * esol_v  = new double [v_nLocBas];
  double * esol_w  = new double [v_nLocBas];

  IQuadPts * quad_vis = new QuadPts_vis_quad4();

  quad_vis -> print_info();

  IQuadPts * quad_gau = new QuadPts_Gauss_Quad( quad_vis->get_num_quadPts_x(), quad_vis->get_num_quadPts_y() );

  quad_gau -> print_info();

  FEAElement * element_quad = new FEAElement_Quad4_3D_der0( quad_vis->get_num_quadPts() );

  // Read the node mappings
  const std::vector<int> analysis_new2old = HDF5_T::read_intVector( "node_mapping_v.h5",
      "/", "new_2_old" );

  double * Rx = new double [v_nLocBas];
  double * Ry = new double [v_nLocBas];
  double * Rz = new double [v_nLocBas];

  // Container for TAWSS & OSI
  std::vector<double> tawss( nFunc, 0.0 ); 
  std::vector<double> osi( nFunc, 0.0 ); 
  std::vector< Vector_3 > osi_top( nFunc, Vector_3(0.0, 0.0, 0.0) );

  const double inv_T = 1.0 / ( static_cast<double>((time_end - time_start)/time_step) + 1.0 );

  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Read solution files
    std::string disp_sol_name(sol_bname), velo_sol_name(sol_bname);
    disp_sol_name.append("disp_");
    velo_sol_name.append("velo_");
    std::string name_to_write(out_bname);

    std::ostringstream time_idx;
    time_idx << 900000000 + time;
    
    disp_sol_name.append(time_idx.str());
    velo_sol_name.append(time_idx.str());
    name_to_write.append(time_idx.str());  

    SYS_T::commPrint("Read %s and %s, and write %s. \n", disp_sol_name.c_str(),
        velo_sol_name.c_str(), name_to_write.c_str() );

    const std::vector<double> disp_sol = ReadPETSc_Vec( disp_sol_name, analysis_new2old, v_nFunc*dof_v, dof_v );
    const std::vector<double> velo_sol = ReadPETSc_Vec( velo_sol_name, analysis_new2old, v_nFunc*dof_v, dof_v );  

    // Ensure that the coordinates of the control points remain in the undeformed state
    v_ctrlPts = v_ctrlPts_origin;
    ctrlPts = ctrlPts_origin;

    // Use the solution displacement field to update the control points
    for(int ii=0; ii<v_nFunc; ++ii)
    {
      v_ctrlPts[3*ii+0] += disp_sol[3*ii+0];
      v_ctrlPts[3*ii+1] += disp_sol[3*ii+1];
      v_ctrlPts[3*ii+2] += disp_sol[3*ii+2];
    }

    for(int ii=0; ii<nFunc; ++ii)
    {
      // use global_node_idx to extract the correct solution entry for this
      // surface node
      ctrlPts[3*ii+0] += disp_sol[3*global_node_idx[ii]+0];
      ctrlPts[3*ii+1] += disp_sol[3*global_node_idx[ii]+1];
      ctrlPts[3*ii+2] += disp_sol[3*global_node_idx[ii]+2];
    }

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

        esol_u[ii] = velo_sol[ dof_v * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 0 ];
        esol_v[ii] = velo_sol[ dof_v * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 1 ];
        esol_w[ii] = velo_sol[ dof_v * v_vecIEN[ee_vol_id*v_nLocBas + ii] + 2 ];
      }

      // Construct the trilinear hexahedron element
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

      std::vector<double> interior_node_coord(3*4*nElem, 0.0); // 4 means that each wall surface element has 4 interior nodes, 3 means x-y-z

      // Now we have found the interior node's volumetric mesh index, record its
      // spatial xyz coordinate
      for (int ii=0; ii<4; ++ii)
      {
        for (int jj=0; jj<3; ++jj)
        {
          interior_node_coord[3*(4*ee + ii )+jj] = v_ctrlPts[ 3*interior_node[4*ee + ii] + jj ];
        }
      }

      // For each nodal point, calculate the outward normal vector using the
      // bilinear quad element. The quad element's ii-th basis corresponds to the
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

      // Loop over the 4 sampling points on the wall bilinear quadrangle element
      for(int qua=0; qua<nLocBas; ++qua)
      {
        // Obtain the 8 basis function's value at the wall boundary points
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

    // Do the averaging by dividing by the area owned by this node
    for(int ii=0; ii<nFunc; ++ii)
    {
      wss_ave[ii].x() /= node_area[ii];
      wss_ave[ii].y() /= node_area[ii];
      wss_ave[ii].z() /= node_area[ii];
    }

    // write the wall shear stress at this time instance
    write_quad_grid_wss( name_to_write, nFunc, nElem, ctrlPts, vecIEN, wss_ave );

    for(int ii=0; ii<nFunc; ++ii)
    {
      tawss[ii]   += inv_T * wss_ave[ii].norm2();
      osi_top[ii] += inv_T * wss_ave[ii];
    } 

  }// Loop over each time instance

  for(int ii=0; ii<nFunc; ++ii)
  {
    const double mag = osi_top[ii].norm2();

    // We will disallow very small tawss in the osi calculation
    if( std::abs(tawss[ii]) > 1.0e-12 )
      osi[ii] = 0.5 * (1.0 - mag / tawss[ii] );
    else
      osi[ii] = 0.0;
  }

  // write time averaged wss and osi
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
  std::vector<int> surface_id_range(4, -1);

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
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, first))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 3;
    surface_id_range[2] = 2;
    surface_id_range[3] = 1;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, second))
  {
    surface_id_range[0] = 2;
    surface_id_range[1] = 3;
    surface_id_range[2] = 7;
    surface_id_range[3] = 6;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, third))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 4;
    surface_id_range[2] = 7;
    surface_id_range[3] = 3;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, fourth))
  {
    surface_id_range[0] = 0;
    surface_id_range[1] = 1;
    surface_id_range[2] = 5;
    surface_id_range[3] = 4;
  }
  else if(std::equal(interior_id_range, interior_id_range + 4, fifth))
  {
    surface_id_range[0] = 1;
    surface_id_range[1] = 2;
    surface_id_range[2] = 6;
    surface_id_range[3] = 5;
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
  vtkPolyData * grid_w = vtkPolyData::New();

  // generate the quad grid
  HEX_T::gen_quad_grid(grid_w, numpts, numcels, pt, ien_array);

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
  vtkPolyData * grid_w = vtkPolyData::New();

  // generate the quad grid
  HEX_T::gen_quad_grid(grid_w, numpts, numcels, pt, ien_array);

  // write tawss
  VTK_T::add_double_PointData(grid_w, tawss, "TAWSS");

  // write osi
  VTK_T::add_double_PointData(grid_w, osi, "OSI");

  // write vtu
  VTK_T::write_vtkPointSet(filename, grid_w, true);

  grid_w->Delete();
}

// EOF
