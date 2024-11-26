#include "ALocal_RotatedBC.hpp"

ALocal_RotatedBC::ALocal_RotatedBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/rotated_nbc"); 

  std::string groupbase(gname);

  Num_LD = h5r -> read_intScalar( gname.c_str(), "Num_LD" );

  // If this sub-domain of this CPU contains local moving bc points, load the LDN array.
  if( Num_LD > 0 )
  {
    LDN = h5r -> read_intVector( gname.c_str(), "LDN" );

    const std::vector<double> temp_ldn_xyz = h5r->read_doubleVector( gname.c_str(), "LDN_pt_xyz" );

    ASSERT( VEC_T::get_size(temp_ldn_xyz) == Num_LD*3, "Error: ALocal_RotatedBC LDN_pt_xyz format is wrong.\n");

    LDN_pt_xyz = std::vector<Vector_3> (Num_LD, Vector_3{ 0, 0, 0 });
    
    for(int ii {0}; ii < Num_LD; ++ii)
      LDN_pt_xyz[ii] = Vector_3{ temp_ldn_xyz[3 * ii], temp_ldn_xyz[3 * ii + 1], temp_ldn_xyz[3 * ii + 2] };
  }
  else
  {
    LDN.clear();
    LDN_pt_xyz.clear();
  }

  SYS_T::print_fatal_if( Num_LD != static_cast<int>( LDN.size() ), "Error: the LDN vector size does not match with the value of Num_LD.\n" );

  // Cell-related data
  cell_nLocBas   = h5r->read_intScalar( gname.c_str(), "cell_nLocBas" );
  num_local_cell = h5r->read_intScalar( gname.c_str(), "num_local_cell" );
  num_local_node = h5r->read_intScalar( gname.c_str(), "num_local_node" );

  // If this partitioned sub-domain contains moving surface element,
  // load its geometrical info 
  if(num_local_cell > 0)
  {      
    local_cell_ien = h5r->read_intVector( gname.c_str(), "local_cell_ien" );
  }
  else
  {
    local_cell_ien.clear();
  }

  if(num_local_node > 0)
  {  
    const std::vector<double> temp_xyz = h5r->read_doubleVector( gname.c_str(), "local_pt_xyz" );

    ASSERT( VEC_T::get_size(temp_xyz) == num_local_node*3, "Error: ALocal_RotatedBC local_pt_xyz format is wrong.\n");

    local_pt_xyz = std::vector<Vector_3> (num_local_node, Vector_3{ 0, 0, 0 });
      
    for(int ii {0}; ii < num_local_node; ++ii)
      local_pt_xyz[ii] = Vector_3{ temp_xyz[3 * ii], temp_xyz[3 * ii + 1], temp_xyz[3 * ii + 2] };

    local_node_pos = h5r->read_intVector( gname.c_str(), "local_node_pos" );
  }
  else
  {
    local_pt_xyz.clear();
    local_node_pos.clear();
  }

  delete h5r; H5Fclose( file_id );    
}

void ALocal_RotatedBC::get_ctrlPts_xyz( const int &eindex, 
    double * const &ctrl_x, 
    double * const &ctrl_y, 
    double * const &ctrl_z ) const
{
  for(int jj=0; jj<cell_nLocBas; ++jj)
  {
    const int pos = local_cell_ien[ cell_nLocBas*eindex+jj ];
    ctrl_x[jj] = local_pt_xyz[pos].x();
    ctrl_y[jj] = local_pt_xyz[pos].y();
    ctrl_z[jj] = local_pt_xyz[pos].z();
  }
}

void ALocal_RotatedBC::get_SIEN( const int &eindex, int * const &sien ) const
{
  for(int jj=0; jj<cell_nLocBas; ++jj)
  {
    const int pos = local_cell_ien[ cell_nLocBas*eindex+jj ];
    sien[jj] = local_node_pos[pos];
  }
}

std::vector<int> ALocal_RotatedBC::get_SIEN( const int &eindex ) const
{
  std::vector<int> out( cell_nLocBas, 0 );
  for(int jj=0; jj<cell_nLocBas; ++jj)
  {
    const int pos = local_cell_ien[ cell_nLocBas*eindex+jj ];
    out[jj] = local_node_pos[pos];
  }
  return out;
}

// EOF
