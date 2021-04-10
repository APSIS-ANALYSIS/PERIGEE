#include "ALocal_Ring_NodalBC.hpp"

ALocal_Ring_NodalBC::ALocal_Ring_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/ring");
  
  Num_LD = h5r -> read_intScalar( gname.c_str(), "Num_LD" );

  num_caps = h5r -> read_intScalar( gname.c_str(), "num_caps" );

  h5r -> read_intVector( gname.c_str(), "cap_dominant_comp", dominant_comp );

  h5r -> read_doubleVector( gname.c_str(), "cap_out_normal", outnormal );

  cap_id.clear();

  // If this sub-domain contains local ring bc points,
  // load the LDN array and corresponding cap ids.
  if( Num_LD > 0 )
  {
    h5r->read_intVector( gname.c_str(), "LDN", LDN );
    h5r->read_intVector( gname.c_str(), "cap_id", cap_id );
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_Ring_NodalBC::~ALocal_Ring_NodalBC()
{
  VEC_T::clean(LDN);
  VEC_T::clean(cap_id);
  VEC_T::clean(dominant_comp);
  VEC_T::clean(outnormal);
}

// EOF
