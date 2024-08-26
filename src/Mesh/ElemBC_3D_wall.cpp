#include "ElemBC_3D_wall.hpp"

ElemBC_3D_wall::ElemBC_3D_wall( const int &elemtype )
: ElemBC_3D( elemtype )
{
  radius.clear();
  thickness.clear();
  youngsmod.clear();
  springconst.clear();
  dampingconst.clear();
}

ElemBC_3D_wall::ElemBC_3D_wall(
    const std::string &walls_combined,
    const double &uniform_thickness,
    const double &uniform_youngsmod,
    const double &uniform_springconst,
    const double &uniform_dampingconst,
    const int &elemtype )
: ElemBC_3D( {walls_combined}, elemtype )
{
  // num_ebc = 1 per the assumption for wall elem bc
  constexpr int ebc_id = 0;

  radius.resize(    num_node[ebc_id] );
  thickness.resize( num_node[ebc_id] );
  youngsmod.resize( num_node[ebc_id] );

  springconst.resize(  num_node[ebc_id] );
  dampingconst.resize( num_node[ebc_id] );

  // Make sure that the files exist
  SYS_T::file_check( walls_combined );

  for(int ii=0; ii<num_node[ebc_id]; ++ii)
  {
    thickness[ii] = uniform_thickness; 
    youngsmod[ii] = uniform_youngsmod;

    springconst[ii]  = uniform_springconst;
    dampingconst[ii] = uniform_dampingconst;

    radius[ii] = 0.0; // radius is not calculated in this case
  }
 
  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");

  std::cout<<"     thickness h = "       << uniform_thickness << std::endl;
  std::cout<<"     Young's modulus E = " << uniform_youngsmod << std::endl;
  std::cout<<"     spring constant ks = "  << uniform_springconst  << std::endl;
  std::cout<<"     damping constant cs = " << uniform_dampingconst << std::endl;
}

ElemBC_3D_wall::ElemBC_3D_wall(
    const std::string &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const double &springconst_combined,
    const double &dampingconst_combined,
    const int &elemtype )
: ElemBC_3D( {walls_combined}, elemtype )
{
  // num_ebc = 1 per the assumption for wall elem bc
  constexpr int ebc_id = 0;

  radius.resize(    num_node[ebc_id] );
  thickness.resize( num_node[ebc_id] );
  youngsmod.resize( num_node[ebc_id] );

  springconst.resize(  num_node[ebc_id] );
  dampingconst.resize( num_node[ebc_id] );

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

  vtkGenericCell * cell = vtkGenericCell::New();

  for(int ii=0; ii<num_node[ebc_id]; ++ii)
  {
    double pt[3] {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

    double cl_pt[3];
    vtkIdType cellId; int subId; double dist;
    
    locator -> FindClosestPoint(&pt[0], &cl_pt[0], cell, cellId, subId, dist); 

    radius[ii] = Vec3::dist( Vector_3(pt[0], pt[1], pt[2]), Vector_3(cl_pt[0], cl_pt[1], cl_pt[2]) );
    thickness[ii] = radius[ii] * thickness2radius_combined; 

    springconst[ii]  = springconst_combined;
    dampingconst[ii] = dampingconst_combined;

    compute_youngsmod(radius[ii], thickness[ii], youngsmod[ii]);
  }
 
  // clean memory
  locator -> Delete();
  reader  -> Delete();
  cell    -> Delete();

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
  
  std::cout<<"     ElemBC_3D_wall generated from "<<walls_combined<<" and ";
  std::cout<<centerlines_combined<<std::endl;
  std::cout<<"     thickness h ranges in ["<<*std::min_element(thickness.begin(), thickness.end())
    <<" , "<<*std::max_element(thickness.begin(), thickness.end())<<"] \n";
  std::cout<<"     Young's modulus E ranges in ["<<*std::min_element(youngsmod.begin(), youngsmod.end())
    <<" , "<<*std::max_element(youngsmod.begin(), youngsmod.end())<<"] \n";
  std::cout<<"     spring constant ks = "  << springconst_combined  << std::endl;
  std::cout<<"     damping constant cs = " << dampingconst_combined << std::endl;
}

ElemBC_3D_wall::ElemBC_3D_wall(
    const std::string &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const double &springconst_combined,
    const double &dampingconst_combined,
    const std::vector<std::string> &wallsList,
    const std::vector<std::string> &centerlinesList,
    const std::vector<double> &thickness2radiusList,
    const std::vector<double> &springconstList,
    const std::vector<double> &dampingconstList,
    const int &elemtype )
: ElemBC_3D_wall( {walls_combined}, centerlines_combined, thickness2radius_combined,
                  springconst_combined, dampingconst_combined, elemtype )
{
  // Check inputs
  SYS_T::print_fatal_if( centerlinesList.size() != wallsList.size(),
    "Error: ElemBC_3D_wall constructor: wallsList and centerlinesList must be of the same length.\n");
  SYS_T::print_fatal_if( thickness2radiusList.size() != wallsList.size(),
    "Error: ElemBC_3D_wall constructor: wallsList and thickness2radiusList must be of the same length.\n");
  SYS_T::print_fatal_if( springconstList.size() != wallsList.size(),
    "Error: ElemBC_3D_wall constructor: wallsList and springconstList must be of the same length.\n");
  SYS_T::print_fatal_if( dampingconstList.size() != wallsList.size(),
    "Error: ElemBC_3D_wall constructor: wallsList and dampingconstList must be of the same length.\n");

  const int num_srfs = static_cast<int>( wallsList.size() );

  for(int ii=0; ii<num_srfs; ++ii)
  {
    SYS_T::file_check(wallsList[ii]);
    SYS_T::file_check(centerlinesList[ii]);
  }
  
  // num_ebc = 1 per the assumption for wall elem bc
  constexpr int ebc_id = 0;
 
  std::cout << "     ===> Overwriting background wall properties in \n";

  // Loop over the surfaces (subsets of the whole wall surface)
  // to locally modify the wall property

  for(int ii=0; ii<num_srfs; ++ii)
  {
    int numpts, numcels;
    std::vector<double> pt;
    std::vector<int> ien_array;

    VTK_T::read_grid( wallsList[ii], numpts, numcels, pt, ien_array );

    const std::vector<int> global_node_idx = VTK_T::read_int_PointData(wallsList[ii], "GlobalNodeID");

    vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
    reader -> SetFileName( centerlinesList[ii].c_str() );
    reader -> Update();
    vtkPolyData * centerlineData = reader -> GetOutput();

    // Identify the closest point (not necessarily a cell vertex) on the centerline
    vtkCellLocator * locator = vtkCellLocator::New();
    locator -> Initialize();
    locator -> SetDataSet( centerlineData );
    locator -> BuildLocator();

    // Data that will be returned by the FindClosestPoint funciton in the
    // for-loop
    vtkGenericCell * cell = vtkGenericCell::New();

    for(int jj=0; jj<numpts; ++jj)
    {
      // Search for corresponding global node ID in walls_combined
      auto it = std::find(global_node[ebc_id].begin(), global_node[ebc_id].end(), global_node_idx[jj]);

      if( it != global_node[ebc_id].end() )
      {
        const int idx = std::distance(global_node[ebc_id].begin(), it);

        double pp[3] {pt_xyz[ebc_id][3*idx], pt_xyz[ebc_id][3*idx+1], pt_xyz[ebc_id][3*idx+2]};

        double cl_pt[3];
        vtkIdType cellId; int subId; double dist;

        locator -> FindClosestPoint(&pp[0], &cl_pt[0], cell, cellId, subId, dist);

        radius[idx] = Vec3::dist( Vector_3(pp[0], pp[1], pp[2]), Vector_3(cl_pt[0], cl_pt[1], cl_pt[2]) );
        thickness[idx] = radius[idx] * thickness2radiusList[ii];

        springconst[idx]  = springconstList[ii];
        dampingconst[idx] = dampingconstList[ii];

        youngsmod[idx] = 7.0e6;
        //compute_youngsmod(radius[idx], thickness[idx], youngsmod[idx]);
      }
      else
        SYS_T::print_fatal( "Error: ElemBC_3D_wall constructor: wallsList does not contain the same global node IDs as walls_combined.\n" );
    }
  
    // clean memory
    locator -> Delete();
    reader  -> Delete();
    cell    -> Delete();

    std::cout << "          " << wallsList[ii] << '\n';
  } // End of loop for ii-th wall surface

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
  
  std::cout<<"     ElemBC_3D_wall generated from "<<walls_combined<<" and ";
  std::cout<<centerlines_combined<<std::endl;
  std::cout<<"     thickness h ranges in ["<<*std::min_element(thickness.begin(), thickness.end())
    <<", "<<*std::max_element(thickness.begin(), thickness.end())<<"] \n";
  std::cout<<"     Young's modulus E ranges in ["<<*std::min_element(youngsmod.begin(), youngsmod.end())
    <<", "<<*std::max_element(youngsmod.begin(), youngsmod.end())<<"] \n";
  std::cout<<"     spring constant ks ranges in ["<<*std::min_element(springconst.begin(), springconst.end())
    <<", "<<*std::max_element(springconst.begin(), springconst.end())<<"] \n";
  std::cout<<"     damping constant cs ranges in ["<<*std::min_element(dampingconst.begin(), dampingconst.end())
    <<", "<<*std::max_element(dampingconst.begin(), dampingconst.end())<<"] \n";
}

void ElemBC_3D_wall::overwrite_from_vtk(
    const std::string &wallprop_vtk,
    const int &type,
    const std::string &vtk_fieldname )
{
  SYS_T::file_check( wallprop_vtk );

  const std::vector<int> global_node_idx = VTK_T::read_int_PointData( wallprop_vtk, "GlobalNodeID");
  const std::vector<double> wallprop = VTK_T::read_double_PointData( wallprop_vtk, vtk_fieldname );
  
  constexpr int ebc_id = 0;
  
  for( int ii = 0; ii < num_node[ebc_id]; ++ii )
  {
    // Search for corresponding global node ID in wallprop_vtk
    auto it = std::find(global_node_idx.begin(), global_node_idx.end(), global_node[ebc_id][ii]);

    if( it != global_node_idx.end() )
    {
      const int idx = std::distance(global_node_idx.begin(), it);

      if( type == 0 )
        thickness[ii] = wallprop[idx];
      else if( type == 1 )
        youngsmod[ii] = wallprop[idx];
      else if( type == 2 )
        springconst[ii] = wallprop[idx];
      else if( type == 3 )
        dampingconst[ii] = wallprop[idx];
      else
        SYS_T::print_fatal("Error: ElemBC_3D_wall::overwrite_from_vtk function: Unknown wallprop type.\n");
    }
  }

  switch( type )
  {
    case 0:
      std::cout << "     ===> Overwriting thickness according to " << wallprop_vtk << std::endl; 
      std::cout << "          thickness h ranges in ["<<*std::min_element(thickness.begin(), thickness.end())
        <<", "<<*std::max_element(thickness.begin(), thickness.end())<<"] \n";
      break;
    case 1:
      std::cout << "     ===> Overwriting Young's modulus E according to " << wallprop_vtk << std::endl; 
      std::cout << "          Young's modulus E ranges in ["<<*std::min_element(youngsmod.begin(), youngsmod.end())
        <<", "<<*std::max_element(youngsmod.begin(), youngsmod.end())<<"] \n";
      break;
    case 2:
      std::cout << "     ===> Overwriting spring constant ks according to " << wallprop_vtk << std::endl; 
      std::cout << "          spring constant ks ranges in ["<<*std::min_element(springconst.begin(), springconst.end())
        <<", "<<*std::max_element(springconst.begin(), springconst.end())<<"] \n";
      break;
    case 3:
      std::cout << "     ===> Overwriting damping constant cs according to " << wallprop_vtk << std::endl; 
      std::cout << "          damping constant cs ranges in ["<<*std::min_element(dampingconst.begin(), dampingconst.end())
        <<", "<<*std::max_element(dampingconst.begin(), dampingconst.end())<<"] \n";
      break;
  }

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
}

void ElemBC_3D_wall::print_info() const
{
  ElemBC_3D::print_info();

  VEC_T::print( radius,       "wall_radius.txt",       '\n');
  VEC_T::print( thickness,    "wall_thickness.txt",    '\n');
  VEC_T::print( youngsmod,    "wall_youngsmod.txt",    '\n');
  VEC_T::print( springconst,  "wall_springconst.txt",  '\n');
  VEC_T::print( dampingconst, "wall_dampingconst.txt", '\n');
}

void ElemBC_3D_wall::write_vtk( const int &ebc_id, const std::string &filename ) const
{
  if(elem_type == 501)
  {
    vtkPolyData * grid_w = vtkPolyData::New();

    TET_T::gen_triangle_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], sur_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    VTK_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete();  
  }
  else if(elem_type == 502)
  {
    vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

    TET_T::gen_quadratic_triangle_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], sur_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    VTK_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete();
  }
  else if(elem_type == 601)
  {
    vtkPolyData * grid_w = vtkPolyData::New();

    HEX_T::gen_quad_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], sur_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    VTK_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete(); 
  }
  else if(elem_type == 602)
  {
    vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

    HEX_T::gen_quadratic_quad_grid( grid_w, num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], sur_ien[ebc_id] );

    add_wall_data(grid_w, ebc_id);

    VTK_T::write_vtkPointSet(filename, grid_w);
    grid_w->Delete();    
  }
  else
    SYS_T::print_fatal("Error: ElemBC_3D_wall::write_vtk fucntion: unknown element type. \n");
}

void ElemBC_3D_wall::add_wall_data( vtkPointSet * const &grid_w, const int &ebc_id ) const
{
  // Add nodal indices
  VTK_T::add_int_PointData( grid_w, global_node[ebc_id], "GlobalNodeID" );

  // Add cell indices
  VTK_T::add_int_CellData( grid_w, global_cell[ebc_id], "GlobalElementID");

  // Add thickness
  VTK_T::add_double_PointData( grid_w, thickness, "Thickness" );

  // Add Young's modulus
  VTK_T::add_double_PointData( grid_w, youngsmod, "YoungsModulus" );

  // Add radius
  VTK_T::add_double_PointData( grid_w, radius, "Radius" );

  // Add spring constant
  VTK_T::add_double_PointData( grid_w, springconst, "SpringConstant" );

  // Add damping constant
  VTK_T::add_double_PointData( grid_w, dampingconst, "DampingConstant" );
}

void ElemBC_3D_wall::compute_youngsmod( const double &r, const double &th, double &E )
{
  //const double alpha = 13.3, beta = 0.3;
  //const double fluid_density = 1.065;
  //const double rho_alpha2 = fluid_density * alpha * alpha;
  //const double beta_exp   = 2.0 * beta - 1.0;
  // Note that r is in mm in Xiao et al. (2013)
  //E = 1.0e4 * rho_alpha2 / ( th * pow( 2.0 * 0.1 * r, beta_exp ) );
  
  E = 1.15e7;
}

// EOF
