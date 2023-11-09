// ============================================================================
// vis_501_fsi_wss.cpp
//
// WSS visualization for 4-node tet elements.
//
// Date: Nov 7 2023
// ============================================================================
#include "HDF5_Tools.hpp"
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap, const int &vec_size, const int &in_dof );

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
  const int dof_v = 3;

  // visualization sampling pattern over time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

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

  // Each surface triangle element requires an additional node: interior_node
  std::vector<int> interior_node( nElem, 0 );

  // Identify the interior node for surface elements
  for(int ee=0; ee<nElem; ++ee)
  {
    std::vector<int> trn(3, 0);
    int ten[4];

    trn[0] = global_node_idx[ vecIEN[3*ee+0] ];
    trn[1] = global_node_idx[ vecIEN[3*ee+1] ];
    trn[2] = global_node_idx[ vecIEN[3*ee+2] ];

    ten[0] = v_vecIEN[ global_ele_idx[ee]*4+0 ];
    ten[1] = v_vecIEN[ global_ele_idx[ee]*4+1 ];
    ten[2] = v_vecIEN[ global_ele_idx[ee]*4+2 ];
    ten[3] = v_vecIEN[ global_ele_idx[ee]*4+3 ];

    int node_check = 0;
    for(int ii=0; ii<4; ++ii)
    {
      if( !VEC_T::is_invec( trn, ten[ii] ) ) interior_node[ee] = ten[ii];
      else node_check += 1;
    }

    SYS_T::print_fatal_if(node_check!=3, "Error: the associated tet element is incompatible with the triangle element.\n");
  }
 
  // Sampling points and element container
  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  // Read the node mappings
  const std::vector<int> analysis_new2old = HDF5_T::read_intVector( "node_mapping_v.h5",
      "/", "new_2_old" );

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

    // Interior nodes's xyz coordinates
    std::vector<double> interior_node_coord( 3*nElem, 0.0 );

    // Unit outward normal vector
    std::vector< Vector_3 > outnormal( nElem, Vector_3(0.0, 0.0, 0.0) );

    // Triangle element surface area
    std::vector<double> tri_area( nElem, 0.0 );

    // Identify the interior node for surface elements
    for(int ee=0; ee<nElem; ++ee)
    {
      std::vector<int> trn(3, 0);

      trn[0] = global_node_idx[ vecIEN[3*ee+0] ];
      trn[1] = global_node_idx[ vecIEN[3*ee+1] ];
      trn[2] = global_node_idx[ vecIEN[3*ee+2] ];

      // Record the interior node's coordinates
      interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
      interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
      interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];

     // Now decide the outward normal vector
      const Vector_3 vec01( v_ctrlPts[3*trn[1]] - v_ctrlPts[3*trn[0]], 
          v_ctrlPts[3*trn[1]+1] - v_ctrlPts[3*trn[0]+1], 
          v_ctrlPts[3*trn[1]+2] - v_ctrlPts[3*trn[0]+2] );

      const Vector_3 vec02( v_ctrlPts[3*trn[2]] - v_ctrlPts[3*trn[0]],
          v_ctrlPts[3*trn[2]+1] - v_ctrlPts[3*trn[0]+1],
          v_ctrlPts[3*trn[2]+2] - v_ctrlPts[3*trn[0]+2] );

      outnormal[ee] = Vec3::cross_product(vec01, vec02); // out = vec01 x vec02
      
      // return out length and scale itself to have unit length
      tri_area[ee] = 0.5 * outnormal[ee].normalize();

      const Vector_3 vec03( v_ctrlPts[interior_node[ee]*3] - v_ctrlPts[3*trn[0]], 
          v_ctrlPts[interior_node[ee]*3+1] - v_ctrlPts[3*trn[0]+1],
          v_ctrlPts[interior_node[ee]*3+2] - v_ctrlPts[3*trn[0]+2] );

      if( Vec3::dot_product(outnormal[ee], vec03)> 0 ) outnormal[ee] *= -1.0;
    }

    std::vector< Vector_3 > wss_ave( nFunc, Vector_3(0.0, 0.0, 0.0) );
    std::vector<double> node_area( nFunc, 0.0 );

    for(int ee=0; ee<nElem; ++ee)
    {
      const int trn[3] { vecIEN[3*ee+0], vecIEN[3*ee+1], vecIEN[3*ee+2] };
      
      const double ectrl_x[4] { ctrlPts[3*trn[0] + 0], ctrlPts[3*trn[1] + 0], ctrlPts[3*trn[2] + 0], interior_node_coord[3*ee + 0] };

      const double ectrl_y[4] { ctrlPts[3*trn[0] + 1], ctrlPts[3*trn[1] + 1], ctrlPts[3*trn[2] + 1], interior_node_coord[3*ee + 1] };

      const double ectrl_z[4] { ctrlPts[3*trn[0] + 2], ctrlPts[3*trn[1] + 2], ctrlPts[3*trn[2] + 2], interior_node_coord[3*ee + 2] };

      element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

      const double esol_u[4] { velo_sol[ global_node_idx[trn[0]] * dof_v ],
        velo_sol[ global_node_idx[trn[1]] * dof_v ],
        velo_sol[ global_node_idx[trn[2]] * dof_v ],
        velo_sol[ interior_node[ee] * dof_v ] };

      const double esol_v[4] { velo_sol[ global_node_idx[trn[0]] * dof_v + 1 ],
        velo_sol[ global_node_idx[trn[1]] * dof_v + 1 ],
        velo_sol[ global_node_idx[trn[2]] * dof_v + 1 ],
        velo_sol[ interior_node[ee] * dof_v + 1 ] };

      const double esol_w[4] { velo_sol[ global_node_idx[trn[0]] * dof_v + 2 ],
        velo_sol[ global_node_idx[trn[1]] * dof_v + 2 ],
        velo_sol[ global_node_idx[trn[2]] * dof_v + 2 ],
        velo_sol[ interior_node[ee] * dof_v + 2 ] };

      const double nx = outnormal[ee].x();
      const double ny = outnormal[ee].y();
      const double nz = outnormal[ee].z();

      for(int qua=0; qua<3; ++qua)
      {
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

        wss_ave[ trn[qua] ].x() += wss_x * tri_area[ee];
        wss_ave[ trn[qua] ].y() += wss_y * tri_area[ee];
        wss_ave[ trn[qua] ].z() += wss_z * tri_area[ee];

        node_area[ trn[qua] ] += tri_area[ee];
      } // Loop over the three surface points 
    } // Loop over surface elements

    // Do the averaging by dividing by the area owned by this node
    for(int ii=0; ii<nFunc; ++ii)
    {
      wss_ave[ii].x() /= node_area[ii];
      wss_ave[ii].y() /= node_area[ii];
      wss_ave[ii].z() /= node_area[ii];
    }

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

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap, const int &vec_size, const int &in_dof )
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
  SYS_T::print_fatal_if( get_sol_temp_size != vec_size, "The solution size %d is not compatible with the size %d given by partition file! \n", get_sol_temp_size, vec_size);

  std::vector<double> veccopy( vec_size, 0.0 );
  
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii) veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // copy the solution varibles to the correct location
  std::vector<double> sol( vec_size, 0.0 );

  // check the nodemap size
  SYS_T::print_fatal_if( VEC_T::get_size(nodemap) * in_dof != vec_size, "Error: node map size is incompatible with the solution length. \n");

  for(int ii=0; ii<VEC_T::get_size(nodemap); ++ii)
  {
    for(int jj=0; jj<in_dof; ++jj)
      sol[in_dof*nodemap[ii]+jj] = veccopy[in_dof*ii+jj];
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

  TET_T::gen_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  VTK_T::add_Vector3_PointData( grid_w, wss_on_node, "WSS" ); 

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
