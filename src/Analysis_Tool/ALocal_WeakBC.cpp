#include "ALocal_WeakBC.hpp"

ALocal_WeakBC::ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank) : weakbc_type {0}
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/weak");

  weakbc_type = h5r -> read_intScalar( gname.c_str(), "weakBC_type" );

  if (weakbc_type > 0)
  {
    num_weak_boundary = h5r -> read_intScalar( gname.c_str(), "num_weak_boundary" );

    C_bI = h5r -> read_doubleScalar( gname.c_str(), "C_bI" );

    num_sur_ele = h5r -> read_intVector( gname.c_str(), "num_local_cell" ); 

    part_vol_ele_id.resize(num_weak_boundary);

    ele_face_id.resize(num_weak_boundary);

    if(weakbc_type == 2)
    {
      num_sur_node = h5r -> read_intVector( gname.c_str(), "num_local_cell_node" );

      sur_node_id.resize(num_weak_boundary);

      Q.resize(num_weak_boundary);
    }

    const std::string groupbase("weakBCid_");

    for(int ii{0}; ii < num_weak_boundary; ++ii)
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      part_vol_ele_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "part_vol_cell_id" );

      ele_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "cell_face_id" );

      if(weakbc_type == 2)
      {
        sur_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "global_node_id" );

        const std::vector<double> temp_Q = h5r -> read_doubleVector( subgroup_name.c_str(), "node_rotation_matrix" );

        for(int node{0}; node < num_sur_node[ii]; ++node)
          Q[ii][node] = Tensor2_3D ( temp_Q[9 * node + 0], temp_Q[9 * node + 1], temp_Q[9 * node + 2],
                                     temp_Q[9 * node + 3], temp_Q[9 * node + 4], temp_Q[9 * node + 5],
                                     temp_Q[9 * node + 6], temp_Q[9 * node + 7], temp_Q[9 * node + 8] ); 
      }
    }
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_WeakBC::~ALocal_WeakBC()
{
  if(weakbc_type > 0)
  {
    VEC_T::clean( num_sur_ele );
    VEC_T::clean( num_sur_node );

    for(int ii{0}; ii < num_weak_boundary; ++ii)
    {
      VEC_T::clean( part_vol_ele_id[ii] );
      VEC_T::clean( ele_face_id[ii] );

      if(weakbc_type == 2)
      {
        VEC_T::clean( sur_node_id[ii] );
        VEC_T::clean( Q[ii] );
      }
    }
  }
}