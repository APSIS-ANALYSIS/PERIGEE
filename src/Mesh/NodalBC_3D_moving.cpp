#include "NodalBC_3D_moving.hpp"

NodalBC_3D_moving::NodalBC_3D_moving( const std::vector<std::string> &inffileList,
    const int &nFunc,
    const int &elemtype )
: num_nbc( static_cast<int>( inffileList.size() ) ), elem_type( elemtype )
{
  // 1. Clear the container for Dirichlet nodes
  dir_nodes.clear();
  num_dir_nodes = 0;

  dir_nodes_on_surface.resize( num_nbc );
  for(int ii=0; ii<num_nbc; ++ii) dir_nodes_on_surface[ii].clear();

  num_dir_nodes_on_surface.resize( num_nbc );

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  sur_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

  // Loop over each surface with id ii
  for( int ii=0; ii<num_nbc; ++ii )
  {
    SYS_T::file_check( inffileList[ii] );

    if( elemtype == 501 )
      nLocBas[ii] = 3;
    else if( elemtype == 502 )
      nLocBas[ii] = 6;
    else if( elemtype == 601 )
      nLocBas[ii] = 4;
    else if( elemtype == 602 )
      nLocBas[ii] = 9;
    else 
      SYS_T::print_fatal("Error: NodalBC_3D_moving::NodalBC_3D_moving: unknown element type.\n");

    VTK_T::read_grid( inffileList[ii], num_node[ii], num_cell[ii], pt_xyz[ii], sur_ien[ii] );

    global_node[ii] = VTK_T::read_int_PointData(inffileList[ii], "GlobalNodeID");
    global_cell[ii] = VTK_T::read_int_CellData(inffileList[ii], "GlobalElementID");  //Nodalbc doesn't record GlobalElementID

    for(unsigned int jj=0; jj<global_node[ii].size(); ++jj)
    {
      SYS_T::print_fatal_if( global_node[ii][jj]<0, "Error: NodalBC_3D_moving::NodalBC_3D_moving: negative nodal index! \n");

      dir_nodes.push_back( global_node[ii][jj] );
      dir_nodes_on_surface[ii].push_back( global_node[ii][jj] );
    }

    num_dir_nodes_on_surface[ii] = dir_nodes_on_surface[ii].size();

    num_dir_nodes = dir_nodes.size();

    VEC_T::sort_unique_resize(dir_nodes);

    SYS_T::print_fatal_if( num_dir_nodes != dir_nodes.size(), "Error: NodalBC_3D_moving::NodalBC_3D_moving: there are repeated nodes in the moving file list.\n" );

    // Generate ID array
    Create_ID( nFunc );
  }

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_moving specified by\n";
  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::cout<<"     nbc_id = "<<ii<<": "<<inffileList[ii]<<" :\n ";
    std::cout<<"          num_node: "<<num_node[ii]<<", num_cell: "<<num_cell[ii]<<'\n';
  }
}