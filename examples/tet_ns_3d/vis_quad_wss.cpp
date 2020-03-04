// ==================================================================
// vis_quad_wss.cpp
//
// WSS visualization for ten-node tet elements.
//
// Date: March 2nd 2020
// ==================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_vis_tri6.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"

void range_generator( const int &ii, std::vector<int> &surface_id_range );

void ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap );

void ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof,
    std::vector<double> &sol );

int get_tri_local_id( const double * const &coor_x,
    const double * const &coor_y,
    const double * const &coor_z,
    const int &len,
    const double &x, const double &y, const double &z,
    const double &tol = 1.0e-8 );

void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< std::vector<double> > &wss_on_node );

void write_triangle_grid_tawss_osi( const std::string &filename,
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

  // take values 0, 1, 2, or 3, it determines the side of the surface in the
  // tetrahedron
  std::vector<int> interior_node_local_index; 
  interior_node_local_index.resize( nElem );

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
    // spatial xyz coordinate
    interior_node_coord[3*ee+0] = v_ctrlPts[ 3*interior_node[ee] + 0 ];
    interior_node_coord[3*ee+1] = v_ctrlPts[ 3*interior_node[ee] + 1 ];
    interior_node_coord[3*ee+2] = v_ctrlPts[ 3*interior_node[ee] + 2 ];
  }

  // Volumetric element visualization sampling point 
  IQuadPts * quad = new QuadPts_vis_tet10_v2();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet10_v2( quad-> get_num_quadPts() );

  double * v_ectrl_x = new double [v_nLocBas];
  double * v_ectrl_y = new double [v_nLocBas];
  double * v_ectrl_z = new double [v_nLocBas];
  double * esol_u  = new double [v_nLocBas];
  double * esol_v  = new double [v_nLocBas];
  double * esol_w  = new double [v_nLocBas];

  IQuadPts * quad_tri_vis = new QuadPts_vis_tri6();

  quad_tri_vis -> print_info();

  IQuadPts * quad_tri_gau = new QuadPts_Gauss_Triangle( quad_tri_vis->get_num_quadPts() );

  quad_tri_gau -> print_info();

  FEAElement * element_tri = new FEAElement_Triangle6_3D_der0( quad_tri_vis-> get_num_quadPts() );

  // Read the mappings of the nodal indices
  std::vector<int> analysis_new2old;
  analysis_new2old.resize(v_nFunc);
  ReadNodeMapping("node_mapping.h5", "new_2_old", v_nFunc, &analysis_new2old[0] );

  double * Rx = new double [v_nLocBas];
  double * Ry = new double [v_nLocBas];
  double * Rz = new double [v_nLocBas];

  // Container for the solution vector
  std::vector<double> sol;

  // Container for (averaged) WSS
  std::vector< std::vector<double> > wss_ave;
  wss_ave.resize( nFunc );
  for(int ii=0; ii<nFunc; ++ii) wss_ave[ii].resize(3);

  // Container for the element area associated with surface nodes
  std::vector<double> node_area; node_area.resize(nFunc);

  // Container for Time averaged WSS and OSI
  std::vector<double> tawss, osi;
  std::vector< std::vector<double> > osi_top;
  tawss.resize( nFunc ); osi.resize( nFunc ); osi_top.resize( nFunc );

  for(int ii=0; ii<nFunc; ++ii)
  {
    tawss[ii]   = 0.0;
    osi[ii]     = 0.0;

    osi_top[ii].resize(3);
    osi_top[ii][0] = 0.0;
    osi_top[ii][1] = 0.0;
    osi_top[ii][2] = 0.0;
  }

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

    PetscPrintf(PETSC_COMM_WORLD, "Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    // Read the solution vector and renumber them based on the nodal mappings
    ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof, sol );

    // Zero the container for the averaged WSS
    for(int ii=0; ii<nFunc; ++ii)
    {
      wss_ave[ii][0] = 0.0;
      wss_ave[ii][1] = 0.0;
      wss_ave[ii][2] = 0.0;

      node_area[ii] = 0.0;
    }

    for(int ee=0; ee<nElem; ++ee)
    {
      // ee element's volumetric element id
      const int ee_vol_id = global_ele_idx[ ee ];

      // get the tetrahedron's control points and associated solution nodal
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

      // Construct the quadratic tetrahedral element
      element -> buildBasis(quad, v_ectrl_x, v_ectrl_y, v_ectrl_z);

      // Obtain the local indices of nodes on the wall surface
      std::vector<int> id_range;
      range_generator( interior_node_local_index[ee], id_range );

      // Obtain the control point coordinates for this element
      double * ectrl_x = new double [nLocBas];
      double * ectrl_y = new double [nLocBas];
      double * ectrl_z = new double [nLocBas];

      for(int ii=0; ii<nLocBas; ++ii)
      {
        ectrl_x[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 0 ];
        ectrl_y[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 1 ];
        ectrl_z[ii] = v_ctrlPts[ 3*v_vecIEN[v_nLocBas * ee_vol_id + id_range[ii] ] + 2 ];
      }

      // Build a basis based on the visualization sampling point for wall
      // triangle element
      element_tri -> buildBasis(quad_tri_vis, ectrl_x, ectrl_y, ectrl_z);

      std::vector< std::vector<double> > outnormal;
      outnormal.resize( nLocBas );

      // For each nodal point, calculate the outward normal vector using the
      // triangle element. The triangle element's ii-th basis corresponds to the
      // tetrahedral element's id_range[ii]-th basis
      for(int ii=0; ii<nLocBas; ++ii)
      {
        double nx, ny, nz, len;

        element_tri -> get_normal_out( ii, ectrl_x[ii], ectrl_y[ii], ectrl_z[ii],
            interior_node_coord[3*ee+0], interior_node_coord[3*ee+1],
            interior_node_coord[3*ee+2], nx, ny, nz, len );

        // id_range[ii] 's outward normal
        outnormal[ii].resize(3);
        outnormal[ii][0] = nx;
        outnormal[ii][1] = ny;
        outnormal[ii][2] = nz;
      }

      // Now calcualte the element surface area
      element_tri -> buildBasis( quad_tri_gau, ectrl_x, ectrl_y, ectrl_z );

      double tri_area = 0.0;
      for(int qua=0; qua<quad_tri_gau->get_num_quadPts(); ++qua)
        tri_area += element_tri->get_detJac(qua) * quad_tri_gau->get_qw(qua);

      for(int ii=0; ii<nLocBas; ++ii)
      {
        ectrl_x[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 0 ];
        ectrl_y[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 1 ];
        ectrl_z[ii] = ctrlPts[ 3*vecIEN[nLocBas * ee + ii ] + 2 ];
      }

      // Loop over the 6 sampling points on the wall quadratic triangle element
      for(int qua=0; qua<nLocBas; ++qua)
      {
        // Obtain the 10 basis function's value at the wall boundary points
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

        // obtain the tet element's id_range[qua] node's outward normal vector 
        const double nx = outnormal[qua][0];
        const double ny = outnormal[qua][1];
        const double nz = outnormal[qua][2];

        const double ax = 2.0 * ux * nx + (uy + vx) * ny + (uz + wx) * nz;
        const double ay = (vx + uy) * nx + 2.0 * vy * ny + (vz + wy) * nz;
        const double az = (wx + uz) * nx + (wy + vz) * ny + 2.0 * wz * nz;

        const double b = ax * nx + ay * ny + az * nz;

        const double wss_x = fluid_mu * ( ax - b * nx );
        const double wss_y = fluid_mu * ( ay - b * ny );
        const double wss_z = fluid_mu * ( az - b * nz );

        const int tri_local_id = get_tri_local_id( ectrl_x, ectrl_y, ectrl_z,
            nLocBas, v_ectrl_x[ id_range[qua] ], v_ectrl_y[ id_range[qua] ], 
            v_ectrl_z[ id_range[qua] ], 1.0e-8 );
       
        const int tri_global_id = vecIEN[nLocBas * ee + tri_local_id];

        wss_ave[ tri_global_id ][0] += wss_x * tri_area;
        wss_ave[ tri_global_id ][1] += wss_y * tri_area;
        wss_ave[ tri_global_id ][2] += wss_z * tri_area;

        node_area[ tri_global_id ] += tri_area;

      } // loop over the sampling points (on surface)

      delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;

    } // End of loop over element

    for(int ii=0; ii<nFunc; ++ii)
    {
      wss_ave[ii][0] /= node_area[ii];
      wss_ave[ii][1] /= node_area[ii];
      wss_ave[ii][2] /= node_area[ii];
    }

    // Write the wall shear stress at this time instance
    write_triangle_grid_wss( name_to_write, nFunc, nElem, ctrlPts, vecIEN, wss_ave );

    for(int ii=0; ii<nFunc; ++ii)
    {
      tawss[ii] += inv_T * std::sqrt( wss_ave[ii][0] * wss_ave[ii][0]
          + wss_ave[ii][1] * wss_ave[ii][1] + wss_ave[ii][2] * wss_ave[ii][2] );

      osi_top[ii][0] += inv_T * wss_ave[ii][0];
      osi_top[ii][1] += inv_T * wss_ave[ii][1];
      osi_top[ii][2] += inv_T * wss_ave[ii][2];
    }

  } // End of loop over time

  MPI_Barrier(PETSC_COMM_WORLD);

  for(int ii=0; ii<nFunc; ++ii)
  {
    const double mag = std::sqrt( osi_top[ii][0] * osi_top[ii][0] +
        osi_top[ii][1] * osi_top[ii][1] + osi_top[ii][2] * osi_top[ii][2] );

    // We will disallow very small tawss in the osi calculation
    if( std::abs(tawss[ii]) > 1.0e-12 )
      osi[ii] = 0.5 * (1.0 - mag / tawss[ii] );
    else
      osi[ii] = 0.0;
  }

  // write the TAWSS and OSI
  std::string tawss_osi_file("SOL_TAWSS_OSI" );
  write_triangle_grid_tawss_osi( tawss_osi_file, nFunc, nElem, ctrlPts, vecIEN, tawss, osi );

  delete [] v_ectrl_x; delete [] v_ectrl_y; delete [] v_ectrl_z;
  delete [] esol_u; delete [] esol_v; delete [] esol_w;
  delete [] Rx; delete [] Ry; delete [] Rz; 
  delete quad; delete element;
  delete quad_tri_vis; delete quad_tri_gau; delete element_tri;

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
      surface_id_range[4] = 9;
      surface_id_range[5] = 8;
      break;
    case 1:
      surface_id_range[0] = 0;
      surface_id_range[1] = 2;
      surface_id_range[2] = 3;
      surface_id_range[3] = 6;
      surface_id_range[4] = 9;
      surface_id_range[5] = 7;
      break;
    case 2:
      surface_id_range[0] = 0;
      surface_id_range[1] = 1;
      surface_id_range[2] = 3;
      surface_id_range[3] = 4;
      surface_id_range[4] = 8;
      surface_id_range[5] = 7;
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


int get_tri_local_id( const double * const &coor_x,
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

  SYS_T::print_fatal("Error in get_tri_local_id.\n");

  return -1;
}

void write_triangle_grid_wss( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector< std::vector<double> > &wss_on_node )
{
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 6*numcels) SYS_T::print_fatal("Error: ien array size does not match the number of cells. \n");

  if(int(wss_on_node.size()) != numpts) SYS_T::print_fatal("Error: wss_on_node size does not match the number of points. \n");

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

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
    vtkQuadraticTriangle * tr = vtkQuadraticTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[6*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[6*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[6*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[6*ii+3] );
    tr->GetPointIds()->SetId( 4, ien_array[6*ii+4] );
    tr->GetPointIds()->SetId( 5, ien_array[6*ii+5] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  
  grid_w -> SetCells(22, cl);
  
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

  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtu");
  writer -> SetFileName( name_to_write.c_str() );
  writer->SetInputData(grid_w);
  writer->Write();

  writer->Delete();
  grid_w->Delete();
}


void write_triangle_grid_tawss_osi( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<double> &tawss,
    const std::vector<double> &osi )
{
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 6*numcels) SYS_T::print_fatal("Error: ien array size does not match the number of cells. \n");

  if(int(tawss.size()) != numpts) SYS_T::print_fatal("Error: tawss size does not match the number of points. \n");

  if(int(osi.size()) != numpts) SYS_T::print_fatal("Error: osi size does not match the number of points. \n");

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

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
    vtkQuadraticTriangle * tr = vtkQuadraticTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[6*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[6*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[6*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[6*ii+3] );
    tr->GetPointIds()->SetId( 4, ien_array[6*ii+4] );
    tr->GetPointIds()->SetId( 5, ien_array[6*ii+5] );
    
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  
  grid_w -> SetCells(22, cl);
  
  cl->Delete();

  // write tawss
  vtkDoubleArray * ptindex = vtkDoubleArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("TAWSS");
  for(int ii=0; ii<numpts; ++ii)
  {
    ptindex -> InsertComponent(ii, 0, tawss[ii]);
  }
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  // write osi
  vtkDoubleArray * vtkosi = vtkDoubleArray::New();
  vtkosi -> SetNumberOfComponents(1);
  vtkosi -> SetName("OSI");
  for(int ii=0; ii<numpts; ++ii)
  {
    vtkosi -> InsertComponent(ii, 0, osi[ii]);
  }
  grid_w -> GetPointData() -> AddArray( vtkosi );
  vtkosi -> Delete();

  // write vtu
  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtu");
  writer -> SetFileName( name_to_write.c_str() );
  writer->SetInputData(grid_w);
  writer->Write();

  writer->Delete();
  grid_w->Delete();
}

// EOF
