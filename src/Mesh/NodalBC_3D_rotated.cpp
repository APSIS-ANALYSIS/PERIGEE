#include "NodalBC_3D_rotated.hpp"

NodalBC_3D_rotated::NodalBC_3D_rotated( const std::string &inffile,
    const int &nFunc,
    const int &elemtype )
: num_nbc( static_cast<int>( inffileList.size() ) ), elem_type( elemtype )
{
  // 1. Clear the container for Dirichlet nodes
  dir_nodes_on_rotated_surface.clear();
  num_dir_nodes_on_rotated_surface = 0;

  // 2. Analyze the file type and read in the data
  num_node = 0;
  num_cell = 0;
  nLocBas = 0;

  sur_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

  SYS_T::file_check( inffile );

  if( elemtype == 501 )
    nLocBas = 3;
  else if( elemtype == 502 )
    nLocBas = 6;
  else if( elemtype == 601 )
    nLocBas = 4;
  else if( elemtype == 602 )
    nLocBas = 9;
  else 
    SYS_T::print_fatal("Error: NodalBC_3D_rotated::NodalBC_3D_rotated: unknown element type.\n");

  VTK_T::read_grid( inffile, num_node, num_cell, pt_xyz, sur_ien );

  global_node = VTK_T::read_int_PointData(inffile, "GlobalNodeID");
  global_cell = VTK_T::read_int_CellData(inffile, "GlobalElementID");  //Nodalbc doesn't record GlobalElementID

  for(unsigned int jj=0; jj<global_node.size(); ++jj)
  {
    SYS_T::print_fatal_if( global_node[jj]<0, "Error: NodalBC_3D_rotated::NodalBC_3D_rotated: negative nodal index! \n");

    dir_nodes_on_rotated_surface.push_back( global_node[jj] );
  }

  num_dir_nodes_on_rotated_surface = dir_nodes_on_rotated_surface.size();

  VEC_T::sort_unique_resize(dir_nodes_on_rotated_surface);

  SYS_T::print_fatal_if( num_dir_nodes_on_rotated_surface != dir_nodes_on_rotated_surface.size(), "Error: NodalBC_3D_rotated::NodalBC_3D_rotated: there are repeated nodes in the rotated file list.\n" );

  // Generate ID array
  Create_ID( nFunc );


  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_rotated specified by\n";
  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::cout<<"     nbc_id = "<<ii<<": "<<inffileList[ii]<<" :\n ";
    std::cout<<"          num_node: "<<num_node[ii]<<", num_cell: "<<num_cell[ii]<<'\n';
  }
}