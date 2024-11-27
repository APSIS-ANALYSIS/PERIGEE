// ============================================================================
// vis_p1_wss.cpp
//
// This is a visualization driver for Wall Shear Stress, which is
// defined on the wall elements only.
// 
// This routine works only for linear tetrahedral element.
//
// Author: Ju Liu
// Date: Sept 16 2019
// ============================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"

std::vector<int> ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size );

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof );

void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< Vector_3 > &wss_on_node );

void write_triangle_grid_tawss_osi( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<double> &tawss,
    const std::vector<double> &osi );

int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  // visualization sampling pattern over time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  constexpr int dof = 4;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: vis_p1_wss is a serial program! \n");

  // Directly read in the volumetric and wall file from the file
  // that record the preprocessor command lines.
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  const std::string geo_file  = cmd_h5r -> read_string("/", "geo_file");
  const std::string wall_file = cmd_h5r -> read_string("/", "sur_file_wall");
  const std::string elemType_str = cmd_h5r -> read_string("/", "elemType");

  const FEType elemType = FE_T::to_FEType(elemType_str);

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Now read the material properties from the solver cmd h5 file
  prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  cmd_h5r = new HDF5_Reader( prepcmd_file );
  
  const double fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // enforce this code is for linear element only
  SYS_T::print_fatal_if( elemType != FEType::Tet4, "Error: element type should be linear tet element.\n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);  

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");

  // Print the key data on screen
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sol_bname: "<<sol_bname<<endl;
  cout<<" -time_start: "<<time_start<<endl;
  cout<<" -time_step: "<<time_step<<endl;
  cout<<" -time_end: "<<time_end<<endl;
  cout<<"----------------------------------\n";
  cout<<" geo_file: "<<geo_file<<endl;
  cout<<" wall_file: "<<wall_file<<endl;
  cout<<" elemType: "<<elemType_str<<endl;
  cout<<" out_bname: "<<out_bname<<endl;
  cout<<" fl_mu: "<<fluid_mu<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Make sure the files exist on disk
  SYS_T::file_check( geo_file );
  SYS_T::file_check( wall_file );

  // Now read in the volumetric mesh info
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Now read the wall surface mesh info
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN;
  
  VTK_T::read_vtp_grid( wall_file, nFunc, nElem, ctrlPts, vecIEN );
  
  const std::vector<int> global_node_idx = VTK_T::read_int_PointData(wall_file, "GlobalNodeID");
  const std::vector<int> global_ele_idx  = VTK_T::read_int_CellData(wall_file, "GlobalElementID");

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";

  // Each surface element requires an additional node. interior_node
  // stores their global indices
  std::vector<int> interior_node( nElem, -1 );

  std::vector<double> interior_node_coord(3*nElem, 0.0);

  // On each surface element, there is a unit outward normal vector
  std::vector< Vector_3 > outnormal( nElem, Vector_3(0.0, 0.0, 0.0) );

  // On each surface element, we store its area
  std::vector<double> tri_area( nElem, 0.0 );

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> trn { global_node_idx[ vecIEN[3*ee+0] ],
      global_node_idx[ vecIEN[3*ee+1] ], global_node_idx[ vecIEN[3*ee+2] ] };
    
    const int ten[4] { v_vecIEN[ global_ele_idx[ee]*4+0 ],
      v_vecIEN[ global_ele_idx[ee]*4+1 ], 
      v_vecIEN[ global_ele_idx[ee]*4+2 ],
      v_vecIEN[ global_ele_idx[ee]*4+3 ] };

    // Locate the interior node's global node index
    int node_check = 0;
    for(int ii=0; ii<4; ++ii)
    {
      const bool gotnode = VEC_T::is_invec( trn, ten[ii] );

      if(!gotnode) interior_node[ee] = ten[ii];
      else node_check += 1;
    }

    SYS_T::print_fatal_if(node_check != 3, "Error: the associated tet element is incompatible with the triangle element.\n");

    // Record the interior node's coordinates
    interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
    interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
    interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];

    // Now decide the outward normal vector
    const Vector_3 l01( v_ctrlPts[3*trn[1]] - v_ctrlPts[3*trn[0]],
        v_ctrlPts[3*trn[1]+1] - v_ctrlPts[3*trn[0]+1],
        v_ctrlPts[3*trn[1]+2] - v_ctrlPts[3*trn[0]+2] );

    const Vector_3 l02( v_ctrlPts[3*trn[2]] - v_ctrlPts[3*trn[0]],
        v_ctrlPts[3*trn[2]+1] - v_ctrlPts[3*trn[0]+1],
        v_ctrlPts[3*trn[2]+2] - v_ctrlPts[3*trn[0]+2] ); 

    Vector_3 ou = Vec3::cross_product( l01, l02 );

    tri_area[ee] = 0.5 * ou.normalize();

    const Vector_3 inw( v_ctrlPts[interior_node[ee]*3] - v_ctrlPts[3*trn[0]],
        v_ctrlPts[interior_node[ee]*3+1] - v_ctrlPts[3*trn[0]+1],
        v_ctrlPts[interior_node[ee]*3+2] - v_ctrlPts[3*trn[0]+2] );

    const double out_dot_in = Vec3::dot_product( ou, inw );

    // if in case the normal points inside, multiply by -1
    if( out_dot_in > 0 ) ou *= -1.0;

    // record the outward normal to outnormal vector 
    outnormal[ee] = ou;
  }

  // Clean the volumetric data for memory
  VEC_T::clean(v_ctrlPts); VEC_T::clean(v_vecIEN);

  // Sampling points for visualization -- classical vertices of tet sampled
  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  // Element container
  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  // Read the node mappings
  const auto analysis_new2old = ReadNodeMapping("node_mapping.h5", "new_2_old", v_nFunc );

  // Read solutions
  std::ostringstream time_index;

  // Container for TAWSS & OSI
  std::vector<double> tawss( nFunc, 0.0 ); 
  std::vector<double> osi( nFunc, 0.0 ); 
  std::vector< Vector_3 > osi_top( nFunc, Vector_3(0.0, 0.0, 0.0) );

  const double inv_T = 1.0 / ( static_cast<double>((time_end - time_start)/time_step) + 1.0 );

  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index << 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    std::cout<<"Time "<<time<<": Read "<<name_to_read<<" and Write "<<name_to_write<<std::endl;

    // Read in the solution vector and arrange them into the natural numbering
    const auto sol = ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof );

    // Container for WSS averaged value
    std::vector< Vector_3 > wss_ave( nFunc, Vector_3(0.0, 0.0, 0.0) );

    // Container for the area associated with the node
    std::vector<double> node_area( nFunc, 0.0 );

    for(int ee=0; ee<nElem; ++ee)
    {
      // Make sure the interior node is the node 3 
      const int trn[3] { vecIEN[3*ee+0], vecIEN[3*ee+1], vecIEN[3*ee+2] };

      const double ectrl_x[4] { ctrlPts[3*trn[0] + 0], ctrlPts[3*trn[1] + 0],
        ctrlPts[3*trn[2] + 0], interior_node_coord[3*ee + 0] };

      const double ectrl_y[4] { ctrlPts[3*trn[0] + 1], ctrlPts[3*trn[1] + 1],
        ctrlPts[3*trn[2] + 1], interior_node_coord[3*ee + 1] };

      const double ectrl_z[4] { ctrlPts[3*trn[0] + 2], ctrlPts[3*trn[1] + 2],
        ctrlPts[3*trn[2] + 2], interior_node_coord[3*ee + 2] };

      element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z); 

      const double esol_u[4] { sol[ global_node_idx[trn[0]] * dof + 1 ],
        sol[ global_node_idx[trn[1]] * dof + 1 ],
        sol[ global_node_idx[trn[2]] * dof + 1 ],
        sol[ interior_node[ee] * dof + 1 ] };

      const double esol_v[4] { sol[ global_node_idx[trn[0]] * dof + 2 ],
        sol[ global_node_idx[trn[1]] * dof + 2 ],
        sol[ global_node_idx[trn[2]] * dof + 2 ],
        sol[ interior_node[ee] * dof + 2 ] };

      const double esol_w[4] { sol[ global_node_idx[trn[0]] * dof + 3 ],
        sol[ global_node_idx[trn[1]] * dof + 3 ],
        sol[ global_node_idx[trn[2]] * dof + 3 ],
        sol[ interior_node[ee] * dof + 3 ] };

      const double nx = outnormal[ee].x();
      const double ny = outnormal[ee].y();
      const double nz = outnormal[ee].z();

      for(int qua=0; qua<3; ++qua)
      {
        // Because the inherent numbering of the element routine, the
        // 4-th basis function is associated with the interior node.
        // Also, the sampling points are aranged such that the first three
        // node are evaluated first. 
        double Rx[4], Ry[4], Rz[4];
        element -> get_gradR(qua, Rx, Ry, Rz);

        const double ux = esol_u[0] * Rx[0] + esol_u[1] * Rx[1] + esol_u[2] * Rx[2] + esol_u[3] * Rx[3];
        const double vx = esol_v[0] * Rx[0] + esol_v[1] * Rx[1] + esol_v[2] * Rx[2] + esol_v[3] * Rx[3];
        const double wx = esol_w[0] * Rx[0] + esol_w[1] * Rx[1] + esol_w[2] * Rx[2] + esol_w[3] * Rx[3];

        const double uy = esol_u[0] * Ry[0] + esol_u[1] * Ry[1] + esol_u[2] * Ry[2] + esol_u[3] * Ry[3];
        const double vy = esol_v[0] * Ry[0] + esol_v[1] * Ry[1] + esol_v[2] * Ry[2] + esol_v[3] * Ry[3];
        const double wy = esol_w[0] * Ry[0] + esol_w[1] * Ry[1] + esol_w[2] * Ry[2] + esol_w[3] * Ry[3];

        const double uz = esol_u[0] * Rz[0] + esol_u[1] * Rz[1] + esol_u[2] * Rz[2] + esol_u[3] * Rz[3];
        const double vz = esol_v[0] * Rz[0] + esol_v[1] * Rz[1] + esol_v[2] * Rz[2] + esol_v[3] * Rz[3];
        const double wz = esol_w[0] * Rz[0] + esol_w[1] * Rz[1] + esol_w[2] * Rz[2] + esol_w[3] * Rz[3];

        const double ax = 2.0 * ux * nx + (uy + vx) * ny + (uz + wx) * nz;
        const double ay = (vx + uy) * nx + 2.0 * vy * ny + (vz + wy) * nz;
        const double az = (wx + uz) * nx + (wy + vz) * ny + 2.0 * wz * nz;

        const double b = ax * nx + ay * ny + az * nz;

        const double wss_x = fluid_mu * ( ax - b * nx );
        const double wss_y = fluid_mu * ( ay - b * ny );
        const double wss_z = fluid_mu * ( az - b * nz ); 

        // Due to the numbering in the visualization quadrature routine, the
        // first quadrature point located on the vertex associated with basis-0,
        // and so on.
        wss_ave[ trn[qua] ] += tri_area[ee] * Vector_3(wss_x, wss_y, wss_z);

        node_area[ trn[qua] ] += tri_area[ee]; 
      } // Loop over the three surface points  
    } // Loop over surface elements

    // Do the averaging by dividing by the area owned by this node
    for(int ii=0; ii<nFunc; ++ii) wss_ave[ii] *= (1.0/ node_area[ii]);

    // write the wall shear stress at this time instance
    write_triangle_grid_wss( name_to_write, nFunc, nElem, ctrlPts, vecIEN, wss_ave );

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
  write_triangle_grid_tawss_osi( tawss_osi_file, nFunc, nElem, ctrlPts, vecIEN, tawss, osi );

  delete quad; delete element;
  PetscFinalize();
  return EXIT_SUCCESS;
}
// END of MAIN Function


// Read the node mappings
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

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space, H5P_DEFAULT, &out[0] );

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

  for(int ii=0; ii<vec_size; ++ii) veccopy[ii] = array_temp[ii];

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

void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< Vector_3 > &wss_on_node )
{
  vtkPolyData * grid_w = vtkPolyData::New();

  // generate the triangle grid
  TET_T::gen_triangle_grid(grid_w, numpts, numcels, pt, ien_array);

  // write wss
  VTK_T::add_Vector3_PointData(grid_w, wss_on_node, "WSS");

  // write vtp
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}

void write_triangle_grid_tawss_osi( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<double> &tawss,
    const std::vector<double> &osi )
{
  vtkPolyData * grid_w = vtkPolyData::New();

  // generate the triangle grid
  TET_T::gen_triangle_grid(grid_w, numpts, numcels, pt, ien_array);

  // write tawss
  VTK_T::add_double_PointData(grid_w, tawss, "TAWSS");

  // write osi 
  VTK_T::add_double_PointData(grid_w, osi, "OSI");
  
  // write vtp
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}

// EOF
