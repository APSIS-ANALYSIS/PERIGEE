#include "APart_Node_wBubble.hpp"

APart_Node_wBubble::APart_Node_wBubble(
    const std::string &fileBaseName, const int &rank,
    const int &nbub_per_cell )
: APart_Node(fileBaseName, rank),
  num_bubble_per_cell( nbub_per_cell ) 
{
  // ----------------------------------------------------------------
  // Read in the number of local cells
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, rank );
  
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  num_local_cell = h5r->read_intScalar( "/Local_Elem", "nlocalele" );

  delete h5r; H5Fclose( file_id );
  // ----------------------------------------------------------------
  // Read in the geo2phy mapper
  hid_t file_id_g2p = H5Fopen("GlobalPartInfo.h5",H5F_ACC_RDONLY,H5P_DEFAULT);
  
  HDF5_Reader * h5r_g2p = new HDF5_Reader( file_id_g2p );

  std::vector<int> g2p;

  h5r_g2p -> read_intVector("/", "geo2phy", g2p);

  delete h5r_g2p; H5Fclose(file_id_g2p);
  // ----------------------------------------------------------------
  // Use ge2p mapper to update the indices for the nodes in
  // local_to_global, node_ghost, node_loc to give them their physics
  // indices
  for(int ii=0; ii<nlocghonode; ++ii)
    local_to_global[ii] = g2p[ local_to_global[ii] ];

  for(int ii=0; ii<nghostnode; ++ii)
    node_ghost[ii] = g2p[ node_ghost[ii] ];

  for(int ii=0; ii<nlocalnode; ++ii)
    node_loc[ii] = g2p[ node_loc[ii] ];
  // ----------------------------------------------------------------

  // Define the total number of bubble enrichment in this cpu patch
  num_total_bubbles = num_local_cell * num_bubble_per_cell;

  // Generate the bubble nodes' indices
  node_bubble.resize( num_total_bubbles );

  // The last local node idx
  const int offset = node_loc[nlocalnode-1] + 1; 
  
  // Define the global bubble indices as following the offset.
  // Their slots have been reserved in preprocessor.
  for(int ee=0; ee<num_local_cell; ++ee)
  {
    for(int ii=0; ii<num_total_bubbles; ++ii)
      node_bubble[ ii ] = offset + ii;
  }

  // Now we append the node_bubble at the end of local_to_global 
  VEC_T::insert_end( local_to_global, node_bubble );

  // Append the node_bubble at the end of node_loc
  VEC_T::insert_end( node_loc, node_bubble );

  // Update the number of total nodes
  ntotalnode = nlocghonode + num_total_bubbles;

  // Update the nmber of local nodes
  nlocalnode += num_total_bubbles;

  // Update the number of locghonode
  nlocghonode += num_total_bubbles;
}


APart_Node_wBubble::~APart_Node_wBubble()
{
  VEC_T::clean( node_bubble );
}


void APart_Node_wBubble::print_info() const
{
  APart_Node::print_info();
  std::cout<<"\n number bubble per cell: "<<num_bubble_per_cell;
  std::cout<<"\n number local cell: "<<num_local_cell;
  std::cout<<"\n number total bubbles: "<<num_total_bubbles;
  std::cout<<"\n node_bubble: \n";
  VEC_T::print( node_bubble );
}

// EOF
