#include "ALocal_EBC_wIntPts.hpp"

ALocal_EBC_wIntPts::ALocal_EBC_wIntPts( const std::string &fileBaseName,
    const int &cpu_rank, const std::string &gname )
: ALocal_EBC( fileBaseName, cpu_rank, gname )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  local_intpts.resize( num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(gname);
      subgroup_name.append("/ebcid_");
      subgroup_name.append( SYS_T::to_string(ii) );

      local_intpts[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "local_intpt_xyz" );
    }
    else
      local_intpts[ii].clear();
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_EBC_wIntPts::~ALocal_EBC_wIntPts()
{
  VEC_T::clean(local_intpts);
}

// EOF
