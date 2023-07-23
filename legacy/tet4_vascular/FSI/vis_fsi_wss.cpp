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

void ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap );

void ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof,
    std::vector<double> &sol );

void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< std::vector<double> > &wss_on_node );

int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  double fluid_mu = 4.0e-2;

  const int dof = 7;

  int time_index = 0;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionInt("-time_index", time_index);

  std::string out_bname = sol_bname;
  out_bname.append("WSS_");

  // Read in the mesh file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  
  const std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string wall_file = cmd_h5r -> read_string("/", "sur_f_file_wall");

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

  SYS_T::file_check( geo_file.c_str() );
  SYS_T::file_check( wall_file.c_str() );

  // ----------------------------------------------------------------
  // Read in the whole FSI volumetric mesh  
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN, phy_tag;
  std::vector<double> v_ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), v_nFunc, v_nElem, v_ctrlPts, v_vecIEN, phy_tag);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Read the wall surface mesh
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN, global_node_idx, global_ele_idx;

  TET_T::read_vtp_grid( wall_file.c_str(), nFunc, nElem, ctrlPts, vecIEN,
      global_node_idx, global_ele_idx );

  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";
  
  // Read the node mappings
  std::vector<int> analysis_new2old;
  analysis_new2old.resize( v_nFunc );
  ReadNodeMapping("node_mapping.h5", "new_2_old", v_nFunc, &analysis_new2old[0] );

  // Solution file name
  std::string name_to_read(sol_bname);
  std::string name_to_write(out_bname);
  std::ostringstream time_idx;
  time_idx.str("");
  time_idx << 900000000 + time_index;
  name_to_read.append(time_idx.str());
  name_to_write.append(time_idx.str());

  PetscPrintf(PETSC_COMM_WORLD, "Read %s and Write %s \n",
      name_to_read.c_str(), name_to_write.c_str() );

  std::vector<double> sol;

  ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof, sol );
  // ----------------------------------------------------------------
  // Use the solution displacement field to update the control points
  for(int ii=0; ii<v_nFunc; ++ii)
  {
    v_ctrlPts[3*ii+0] += sol[7*ii+0];
    v_ctrlPts[3*ii+1] += sol[7*ii+1];
    v_ctrlPts[3*ii+2] += sol[7*ii+2];
  }

  for(int ii=0; ii<nFunc; ++ii)
  {
    int idx = global_node_idx[ii];
    ctrlPts[3*ii+0] += sol[7*idx+0];
    ctrlPts[3*ii+1] += sol[7*idx+1];
    ctrlPts[3*ii+2] += sol[7*idx+2];
  }

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

  // Identify the interior node for surface elements
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
    
    // Record the interior node's coordinates
    interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
    interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
    interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];
 
   // Now decide the outward normal vector
    outnormal[ee].resize(3);
    const double l01x = v_ctrlPts[3*trn[1]] - v_ctrlPts[3*trn[0]];
    const double l01y = v_ctrlPts[3*trn[1]+1] - v_ctrlPts[3*trn[0]+1];
    const double l01z = v_ctrlPts[3*trn[1]+2] - v_ctrlPts[3*trn[0]+2];

    const double l02x = v_ctrlPts[3*trn[2]] - v_ctrlPts[3*trn[0]];
    const double l02y = v_ctrlPts[3*trn[2]+1] - v_ctrlPts[3*trn[0]+1];
    const double l02z = v_ctrlPts[3*trn[2]+2] - v_ctrlPts[3*trn[0]+2];

    double oux, ouy, ouz;

    MATH_T::cross3d( l01x, l01y, l01z, l02x, l02y, l02z, oux, ouy, ouz );

    tri_area[ee] = 0.5 * MATH_T::normalize3d( oux, ouy, ouz );

    const double inwx = v_ctrlPts[interior_node[ee]*3]   - v_ctrlPts[3*trn[0]];
    const double inwy = v_ctrlPts[interior_node[ee]*3+1] - v_ctrlPts[3*trn[0]+1];
    const double inwz = v_ctrlPts[interior_node[ee]*3+2] - v_ctrlPts[3*trn[0]+2];

    const double out_dot_in = MATH_T::dot3d(oux, ouy, ouz, inwx, inwy, inwz); 
    
    if(out_dot_in > 0)
    {
      oux *= -1.0;
      ouy *= -1.0;
      ouz *= -1.0;
    }

    outnormal[ee][0] = oux;
    outnormal[ee][1] = ouy;
    outnormal[ee][2] = ouz;
  }

  // Clean the volumetric data to save memory
  VEC_T::clean(v_ctrlPts); VEC_T::clean(v_vecIEN);

  // Sampling points and element container
  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );
  
  double ectrl_x[4], ectrl_y[4], ectrl_z[4];
  double esol_u[4], esol_v[4], esol_w[4];
  double Rx[4], Ry[4], Rz[4];

  std::vector< std::vector<double> > wss_ave;
  wss_ave.resize( nFunc );
  for(int ii=0; ii<nFunc; ++ii) wss_ave[ii].resize(3);

  std::vector<double> node_area; node_area.resize(nFunc);

  for(int ii=0; ii<nFunc; ++ii)
  {
    wss_ave[ii][0] = 0.0;
    wss_ave[ii][1] = 0.0;
    wss_ave[ii][2] = 0.0;

    node_area[ii] = 0.0;
  }
  
  for(int ee=0; ee<nElem; ++ee)
  {
    double trn[3];
    trn[0] = vecIEN[3*ee+0];
    trn[1] = vecIEN[3*ee+1];
    trn[2] = vecIEN[3*ee+2];

    ectrl_x[0] = ctrlPts[3*trn[0] + 0];
    ectrl_x[1] = ctrlPts[3*trn[1] + 0];
    ectrl_x[2] = ctrlPts[3*trn[2] + 0];
    ectrl_x[3] = interior_node_coord[3*ee + 0];

    ectrl_y[0] = ctrlPts[3*trn[0] + 1];
    ectrl_y[1] = ctrlPts[3*trn[1] + 1];
    ectrl_y[2] = ctrlPts[3*trn[2] + 1];
    ectrl_y[3] = interior_node_coord[3*ee + 1];

    ectrl_z[0] = ctrlPts[3*trn[0] + 2];
    ectrl_z[1] = ctrlPts[3*trn[1] + 2];
    ectrl_z[2] = ctrlPts[3*trn[2] + 2];
    ectrl_z[3] = interior_node_coord[3*ee + 2];

    element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

    esol_u[0] = sol[ global_node_idx[trn[0]] * dof + 4 ];
    esol_u[1] = sol[ global_node_idx[trn[1]] * dof + 4 ];
    esol_u[2] = sol[ global_node_idx[trn[2]] * dof + 4 ];
    esol_u[3] = sol[ interior_node[ee] * dof + 4 ];

    esol_v[0] = sol[ global_node_idx[trn[0]] * dof + 5 ];
    esol_v[1] = sol[ global_node_idx[trn[1]] * dof + 5 ];
    esol_v[2] = sol[ global_node_idx[trn[2]] * dof + 5 ];
    esol_v[3] = sol[ interior_node[ee] * dof + 5 ];

    esol_w[0] = sol[ global_node_idx[trn[0]] * dof + 6 ];
    esol_w[1] = sol[ global_node_idx[trn[1]] * dof + 6 ];
    esol_w[2] = sol[ global_node_idx[trn[2]] * dof + 6 ];
    esol_w[3] = sol[ interior_node[ee] * dof + 6 ];

    const double nx = outnormal[ee][0];
    const double ny = outnormal[ee][1];
    const double nz = outnormal[ee][2];

    for(int qua=0; qua<3; ++qua)
    {
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

      wss_ave[ trn[qua] ][0] += wss_x * tri_area[ee];
      wss_ave[ trn[qua] ][1] += wss_y * tri_area[ee];
      wss_ave[ trn[qua] ][2] += wss_z * tri_area[ee];

      node_area[ trn[qua] ] += tri_area[ee];
    }
  }

  for(int ii=0; ii<nFunc; ++ii)
  {
    wss_ave[ii][0] /= node_area[ii];
    wss_ave[ii][1] /= node_area[ii];
    wss_ave[ii][2] /= node_area[ii];
  }

  write_triangle_grid_wss( name_to_write, nFunc, nElem, ctrlPts, vecIEN, wss_ave );

  MPI_Barrier(PETSC_COMM_WORLD);
  delete quad; delete element;
  PetscFinalize();
  return EXIT_SUCCESS;
}


// Read the node mappings
void ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap )
{
  hid_t file_id = H5Fopen(node_mapping_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT);

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0];

  if( int(dSize) != node_size )
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, nodemap );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);
}

void ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof,
    std::vector<double> &sol )
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
    PetscPrintf(PETSC_COMM_SELF,
        "The solution size %d is not compatible with the size %d given by partition file! \n",
        get_sol_temp_size, vec_size);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<double> veccopy;
  veccopy.resize(vec_size);
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii)
    veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // copy the solution varibles to the correct location
  sol.clear();
  sol.resize(vec_size);

  // check the nodemap size
  if( (int)nodemap.size() * in_dof != vec_size ) SYS_T::print_fatal("Error: node map size is incompatible with the solution length. \n");

  for(unsigned int ii=0; ii<nodemap.size(); ++ii)
  {
    const int index = nodemap[ii];
    for(int jj=0; jj<in_dof; ++jj)
      sol[in_dof*index+jj] = veccopy[in_dof*ii+jj];
  }
}


void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< std::vector<double> > &wss_on_node )
{
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 3*numcels) SYS_T::print_fatal("Error: ien array size does not match the number of cells. \n");

  if(int(wss_on_node.size()) != numpts) SYS_T::print_fatal("Error: wss_on_node size does not match the number of points. \n");

  vtkPolyData * grid_w = vtkPolyData::New();

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. cell data
  vtkCellArray * cl = vtkCellArray::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkTriangle * tr = vtkTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[3*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[3*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[3*ii+2] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  grid_w->SetPolys(cl);
  cl->Delete();

  // write wss
  vtkDoubleArray * ptindex = vtkDoubleArray::New();
  ptindex -> SetNumberOfComponents(3);
  ptindex -> SetName("WSS");
  for(int ii=0; ii<numpts; ++ii)
  {
    ptindex -> InsertComponent(ii, 0, wss_on_node[ii][0]);
    ptindex -> InsertComponent(ii, 1, wss_on_node[ii][1]);
    ptindex -> InsertComponent(ii, 2, wss_on_node[ii][2]);
  }
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtp");
  writer -> SetFileName( name_to_write.c_str() );
  writer->SetInputData(grid_w);
  writer->Write();

  writer->Delete();
  grid_w->Delete();
}

// EOF