#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::string &walls_combined,
    const double &uniform_thickness,
    const double &uniform_youngsmod,
    const int &elemtype,
    const double &in_fluid_density )
: ElemBC_3D_tet( walls_combined, elemtype ),
  elemType( elemtype ),
  fluid_density( in_fluid_density )
{
  // num_ebc = 1 per the assumption for wall elem bc
  const int ebc_id = 0;

  radius.resize(    num_node[ebc_id] );
  thickness.resize( num_node[ebc_id] );
  youngsmod.resize( num_node[ebc_id] );

  // Make sure that the files exist
  SYS_T::file_check( walls_combined );

  for(int ii=0; ii<num_node[ebc_id]; ++ii)
  {
    thickness[ii] = uniform_thickness; 

    youngsmod[ii] = uniform_youngsmod;

    radius[ii] = 0.0; // radius is not calculated in this case
  }
 
  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
}

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::string &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const int &elemtype,
    const double &in_fluid_density )
: ElemBC_3D_tet( walls_combined, elemtype ),
  elemType( elemtype ),
  fluid_density( in_fluid_density )
{
  // num_ebc = 1 per the assumption for wall elem bc
  const int ebc_id = 0;

  radius.resize(    num_node[ebc_id] );
  thickness.resize( num_node[ebc_id] );
  youngsmod.resize( num_node[ebc_id] );

  // Make sure that the files exist
  SYS_T::file_check( walls_combined );
  SYS_T::file_check( centerlines_combined );

  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( centerlines_combined.c_str() );
  reader -> Update();
  vtkPolyData * centerlineData = reader -> GetOutput();

  // Identify the closest point (not necessarily a cell vertex) on the centerline
  vtkCellLocator * locator = vtkCellLocator::New();
  locator -> Initialize();
  locator -> SetDataSet( centerlineData );
  locator -> BuildLocator();

  double cl_pt[3];
  vtkGenericCell * cell = vtkGenericCell::New();
  vtkIdType cellId; int subId; double dist;

  for(int ii=0; ii<num_node[ebc_id]; ++ii)
  {
    double pt[3] = {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

    locator -> FindClosestPoint(&pt[0], &cl_pt[0], cell, cellId, subId, dist); 

    radius[ii] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);
  
    thickness[ii] = radius[ii] * thickness2radius_combined; 

    compute_youngsmod(radius[ii], thickness[ii], youngsmod[ii]);
  }
 
  // clean memory
  locator -> Delete();
  reader  -> Delete();
  cell    -> Delete();

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
}


ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::string &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const std::vector<std::string> &wallsList,
    const std::vector<std::string> &centerlinesList,
    const std::vector<double> &thickness2radiusList,
    const int &elemtype,
    const double &in_fluid_density )
: ElemBC_3D_tet_wall( walls_combined, centerlines_combined, thickness2radius_combined,
                      elemtype, in_fluid_density)
{
  // Check inputs
  SYS_T::print_fatal_if( centerlinesList.size() != wallsList.size(),
    "ERROR: wallsList and centerlinesList must be of the same length.\n");
  SYS_T::print_fatal_if( thickness2radiusList.size() != wallsList.size(),
    "ERROR: wallsList and thickness2radiusList must be of the same length.\n");

  const int num_srfs = static_cast<int>( wallsList.size() );

  for(int ii=0; ii<num_srfs; ++ii)
  {
    SYS_T::file_check(wallsList[ii]);
    SYS_T::file_check(centerlinesList[ii]);
  }
  
  // num_ebc = 1 per the assumption for wall elem bc
  const int ebc_id = 0;
 
  std::cout << "     ===> Overwriting background wall properties in \n";

  // Loop over the surfaces (subsets of the whole wall surface)
  // to locally modify the wall property
  int numpts, numcels;
  std::vector<double> pt;
  std::vector<int> ien_array, global_node_idx, global_elem_idx;
  for(int ii=0; ii<num_srfs; ++ii)
  {
    if(elemtype == 501)
    {
      TET_T::read_vtp_grid( wallsList[ii], numpts, numcels,
            pt, ien_array, global_node_idx, global_elem_idx ); 
    }
    else
    {
      TET_T::read_vtu_grid( wallsList[ii], numpts, numcels,
            pt, ien_array, global_node_idx, global_elem_idx );
    }

    vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
    reader -> SetFileName( centerlinesList[ii].c_str() );
    reader -> Update();
    vtkPolyData * centerlineData = reader -> GetOutput();

    // Identify the closest point (not necessarily a cell vertex) on the centerline
    vtkCellLocator * locator = vtkCellLocator::New();
    locator -> Initialize();
    locator -> SetDataSet( centerlineData );
    locator -> BuildLocator();

    double cl_pt[3];
    vtkGenericCell * cell = vtkGenericCell::New();
    vtkIdType cellId; int subId; double dist;

    for(int jj=0; jj<numpts; ++jj)
    {
      // Search for corresponding global node ID in walls_combined
      auto it = std::find(global_node[ebc_id].begin(), global_node[ebc_id].end(), global_node_idx[jj]);

      if( it != global_node[ebc_id].end() )
      {
        const int idx = std::distance(global_node[ebc_id].begin(), it);

        double pt[3] = {pt_xyz[ebc_id][3*idx], pt_xyz[ebc_id][3*idx+1], pt_xyz[ebc_id][3*idx+2]};

        locator -> FindClosestPoint(&pt[0], &cl_pt[0], cell, cellId, subId, dist);

        radius[idx] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);

        thickness[idx] = radius[idx] * thickness2radiusList[ii];

        compute_youngsmod(radius[idx], thickness[idx], youngsmod[idx]);
      }
      else
        SYS_T::print_fatal( "Error: wallsList does not contain the same global node IDs as walls_combined.\n" );
    }
  
    // clean memory
    locator -> Delete();
    reader  -> Delete();
    cell    -> Delete();

    std::cout << "          " << wallsList[ii] << '\n';
  } // End of loop for ii-th wall surface

  VEC_T::clean( pt ); VEC_T::clean( ien_array ); 
  VEC_T::clean( global_node_idx ); VEC_T::clean( global_elem_idx );

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius    );
  VEC_T::clean( thickness );
  VEC_T::clean( youngsmod );
}


void ElemBC_3D_tet_wall::print_info() const
{
  ElemBC_3D_tet::print_info();

  VEC_T::print( radius,    "wall_radius.txt",    '\n');
  VEC_T::print( thickness, "wall_thickness.txt", '\n');
  VEC_T::print( youngsmod, "wall_youngsmod.txt", '\n');
}


void ElemBC_3D_tet_wall::write_vtk( const int &ebc_id,
    const std::string &filename ) const
{
  if(elemType == 501)
  {
    vtkPolyData * grid_w = vtkPolyData::New();

    TET_T::gen_triangle_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], tri_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    TET_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete();  
  }
  else
  {
    vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

    TET_T::gen_quadratic_triangle_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], tri_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    TET_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete();
  }
}


void ElemBC_3D_tet_wall::add_wall_data( vtkPointSet * const &grid_w, const int &ebc_id ) const
{
  // Add nodal indices
  TET_T::add_int_PointData( grid_w, global_node[ebc_id], "GlobalNodeID" );

  // Add cell indices
  TET_T::add_int_CellData( grid_w, global_cell[ebc_id], "GlobalElementID");

  // Add thickness
  TET_T::add_double_PointData( grid_w, thickness, "Thickness" );

  // Add Young's modulus
  TET_T::add_double_PointData( grid_w, youngsmod, "YoungsModulus" );

  // Add radius
  TET_T::add_double_PointData( grid_w, radius, "Radius" );
}


void ElemBC_3D_tet_wall::compute_youngsmod( const double &r, const double &th, double &E )
{
  const double alpha = 13.3, beta = 0.3;
  const double rho_alpha2 = fluid_density * alpha * alpha;
  const double beta_exp   = 2.0 * beta - 1.0;

  // Note that r is in mm in Xiao et al. (2013)
  E = 1.0e4 * rho_alpha2 / ( th * pow( 2.0 * 0.1 * r, beta_exp ) );
}

// EOF
