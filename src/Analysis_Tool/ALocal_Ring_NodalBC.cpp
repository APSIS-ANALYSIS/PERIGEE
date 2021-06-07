#include "ALocal_Ring_NodalBC.hpp"

ALocal_Ring_NodalBC::ALocal_Ring_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/ring");
  
  ringbc_type = h5r -> read_intScalar( gname.c_str(), "ring_bc_type" );

  Num_LD = h5r -> read_intScalar( gname.c_str(), "Num_LD" );

  num_caps = h5r -> read_intScalar( gname.c_str(), "num_caps" );

  h5r -> read_intVector( gname.c_str(), "cap_dominant_n_comp", dominant_n_comp );

  const std::vector<double> outnormal_vec = h5r -> read_doubleVector( gname.c_str(), "cap_out_normal" );

  outnormal.resize(num_caps);
  for(int ii=0; ii<num_caps; ++ii)
  {
    outnormal[ii](0) = outnormal_vec[3*ii+0];
    outnormal[ii](1) = outnormal_vec[3*ii+1];
    outnormal[ii](2) = outnormal_vec[3*ii+2];
  }

  // If this sub-domain contains local ring nodes,
  // load the LDN array, the corresponding cap ids, unit tangential vectors,
  // and the dominant tangential components 
  if( Num_LD > 0 )
  {
    LDN = h5r->read_intVector( gname.c_str(), "LDN" );
    local_cap_id = h5r->read_intVector( gname.c_str(), "local_cap_id" );
    local_dominant_t_comp = h5r->read_intVector( gname.c_str(), "local_dominant_t_comp" );

    const std::vector<double> tangential_vec = h5r->read_doubleVector( gname.c_str(), "local_tangential" );

    local_tangential.resize(Num_LD);
    for(int ii=0; ii<Num_LD; ++ii)
    { 
      local_tangential[ii](0) = tangential_vec[3*ii+0];
      local_tangential[ii](1) = tangential_vec[3*ii+1];
      local_tangential[ii](2) = tangential_vec[3*ii+2];
    }
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_Ring_NodalBC::~ALocal_Ring_NodalBC()
{
  VEC_T::clean( LDN                   );
  VEC_T::clean( local_cap_id          );
  VEC_T::clean( dominant_n_comp       );
  VEC_T::clean( local_dominant_t_comp );
  VEC_T::clean( outnormal             );
  VEC_T::clean( local_tangential      );
}

// EOF
