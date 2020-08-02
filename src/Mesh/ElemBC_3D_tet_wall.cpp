#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const std::vector<std::string> &centerlineList,
    const std::vector<double> &thickness_to_radius,
    const std::vector<double> &youngsmod_alpha,
    const std::vector<double> &youngsmod_beta,
    const double &fluid_density,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  // Check inputs
  SYS_T::print_fatal_if( centerlineList.size() != vtkfileList.size(),
      "Error: centerlineList length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( thickness_to_radius.size() != vtkfileList.size(),
      "Error: thickness_to_radius length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( youngsmod_alpha.size() != vtkfileList.size(),
      "Error: youngsmod_alpha length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( youngsmod_beta.size() != vtkfileList.size(),
      "Error: youngsmod_beta length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( elemtype != 501,
      "Error: unsupported element type.\n");

  radius.resize(num_ebc);
  thickness.resize(num_ebc);
  youngsmod.resize(num_ebc);
  for(int ii=0; ii<num_ebc; ++ii) 
  {
    radius[ii].resize( num_node[ii] );
    thickness[ii].resize( num_node[ii] );
    youngsmod[ii].resize( num_node[ii] );
  }

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const double rho_alpha2 = fluid_density * youngsmod_alpha[ebc_id] * youngsmod_alpha[ebc_id];
    const double beta_exp   = 2.0 * youngsmod_beta[ebc_id] - 1.0; 

    vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
    reader -> SetFileName( centerlineList[ebc_id].c_str() );
    reader -> Update();
    vtkPolyData * centerlineData = reader -> GetOutput();

    vtkPointLocator * locator = vtkPointLocator::New();
    locator -> Initialize();
    locator -> SetDataSet( centerlineData );
    locator -> BuildLocator();

    for(int ii=0; ii<num_node[ebc_id]; ++ii)
    {
      const double pt[3] = {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

      const int closest_id = locator -> FindClosestPoint(&pt[0]);

      const double * cl_pt = centerlineData -> GetPoints() -> GetPoint(closest_id);

      radius[ebc_id][ii] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);
   
      thickness[ebc_id][ii] = radius[ebc_id][ii] * thickness_to_radius[ebc_id]; 

      youngsmod[ebc_id][ii] = rho_alpha2 / ( thickness[ebc_id][ii] * pow( 2.0*radius[ebc_id][ii], beta_exp ) );
    }

    // clean memory
    locator -> Delete();
    reader -> Delete();
  }

  // Write out vtp's with wall properties
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
    write_wall_prop(ebc_id, "varwallprop_" + SYS_T::to_string(ebc_id));
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  for(int ii=0; ii<num_ebc; ++ii)
  {
    VEC_T::clean( radius[ii]    );
    VEC_T::clean( thickness[ii] );
    VEC_T::clean( youngsmod[ii] );
  }

  VEC_T::clean( radius    );
  VEC_T::clean( thickness );
  VEC_T::clean( youngsmod );
}


void ElemBC_3D_tet_wall::print_info() const
{
  ElemBC_3D_tet::print_info();

  for(int face=0; face<num_ebc; ++face)
  {
    VEC_T::print( radius[face],    "wall_id_" + SYS_T::to_string(face) + "_radius.txt",    '\n');
    VEC_T::print( thickness[face], "wall_id_" + SYS_T::to_string(face) + "_thickness.txt", '\n');
    VEC_T::print( youngsmod[face], "wall_id_" + SYS_T::to_string(face) + "_youngsmod.txt", '\n');
  }
}


void ElemBC_3D_tet_wall::write_vtk( const int &ebc_id, const std::string &filename ) const
{
  vtkPolyData * grid_w = vtkPolyData::New();
  
  TET_T::gen_triangle_grid( grid_w, 
      num_node[ebc_id], num_cell[ebc_id],
      pt_xyz[ebc_id], tri_ien[ebc_id] );    

  // Add nodal indices
  TET_T::add_int_PointData( grid_w, global_node[ebc_id], "GlobalNodeID" );

  // Add cell indices
  TET_T::add_int_CellData( grid_w, global_cell[ebc_id], "GlobalElementID");

  // Add thickness
  TET_T::add_double_PointData( grid_w, thickness[ebc_id], "Thickness" );

  // Add Young's modulus
  TET_T::add_double_PointData( grid_w, youngsmod[ebc_id], "YoungsModulus" );

  // write vtp
  TET_T::write_vtkPointSet(filename, grid_w);

  // Clean memory
  grid_w->Delete();
}

// EOF
