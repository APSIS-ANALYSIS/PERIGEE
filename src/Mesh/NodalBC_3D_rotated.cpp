#include "NodalBC_3D_rotated.hpp"

NodalBC_3D_rotated::NodalBC_3D_rotated( 
    const std::string &rotated_file,
    const std::string &fixed_file,
    const int &nFunc, const FEType &elemtype )
: elem_type( elemtype )
{
  // Prepare the numbers that need to be shifted
  const int fixed_nFunc = VTK_T::read_num_pt(fixed_file);
  const int fixed_nElem = VTK_T::read_num_cl(fixed_file);

  // Clear the container for Dirichlet nodes
  dir_nodes_on_rotated_surface.clear();
  dir_nodes.clear();
  num_dir_nodes_on_rotated_surface = 0;

  // Analyze the file type and read in the data
  num_node = 0;
  num_cell = 0;
  nLocBas = 0;

  SYS_T::file_check( rotated_file );
  SYS_T::file_check( fixed_file );

  if( elem_type == FEType::Tet4 )
    nLocBas = 3;
  else if( elem_type == FEType::Tet10 )
    nLocBas = 6;
  else if( elem_type == FEType::Hex8 )
    nLocBas = 4;
  else if( elem_type == FEType::Hex27 )
    nLocBas = 9;
  else 
    SYS_T::print_fatal("Error: NodalBC_3D_rotated::NodalBC_3D_rotated: unknown element type.\n");

  VTK_T::read_grid( rotated_file, num_node, num_cell, pt_xyz, sur_ien );

  global_node = VTK_T::read_int_PointData(rotated_file, "GlobalNodeID");
  global_cell = VTK_T::read_int_CellData(rotated_file, "GlobalElementID");

  for(int &nodeid : global_node)
    nodeid += fixed_nFunc;

  for(int &cellid : global_cell)
    cellid += fixed_nElem;

  for(unsigned int jj=0; jj<global_node.size(); ++jj)
  {
    SYS_T::print_fatal_if( global_node[jj]<0, "Error: NodalBC_3D_rotated::NodalBC_3D_rotated: negative nodal index! \n");

    dir_nodes_on_rotated_surface.push_back( global_node[jj] );

    dir_nodes.push_back( global_node[jj] );
  }

  num_dir_nodes_on_rotated_surface = dir_nodes_on_rotated_surface.size();

  num_dir_nodes = dir_nodes.size();
 
  VEC_T::sort_unique_resize(dir_nodes_on_rotated_surface);

  VEC_T::sort_unique_resize(dir_nodes);

  SYS_T::print_fatal_if( num_dir_nodes_on_rotated_surface != dir_nodes_on_rotated_surface.size(), "Error: NodalBC_3D_rotated::NodalBC_3D_rotated: there are repeated nodes in the rotated file list.\n" );

  // Generate ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_rotated specified by\n";
  std::cout<<"     "<<rotated_file<<" :\n ";
  std::cout<<"          num_node: "<<num_node<<", num_cell: "<<num_cell<<'\n';
}
