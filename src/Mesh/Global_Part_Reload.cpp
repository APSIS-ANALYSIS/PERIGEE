#include "Global_Part_Reload.hpp"

Global_Part_Reload::Global_Part_Reload( const std::string &element_part_name,
    const std::string &node_part_name )
{
  // --------------------------------------------------------------------------
  const std::string efName = element_part_name + ".h5";

  hid_t efile_id = H5Fopen( efName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * eh5r = new HDF5_Reader( efile_id );
  
  int temp = eh5r -> read_intScalar("/", "isMETIS");
  if( temp == 1 ) isMETIS = true;
  else isMETIS = false;

  temp = eh5r -> read_intScalar("/", "part_isdual");
  if(temp == 1) isDual = true;
  else isDual = false;

  temp = eh5r -> read_intScalar("/", "isSerial");
  if(temp == 1) isSerial = true;
  else isSerial = false;

  dual_edge_ncommon = eh5r -> read_intScalar("/", "in_ncommon");

  epart = eh5r -> read_intVector("/", "part");

  delete eh5r; H5Fclose( efile_id );

  // --------------------------------------------------------------------------
  const std::string nfName = node_part_name + ".h5";

  hid_t nfile_id = H5Fopen( nfName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * nh5r = new HDF5_Reader( nfile_id );

  npart = nh5r -> read_intVector("/", "part");

  delete nh5r; H5Fclose( nfile_id );
  // --------------------------------------------------------------------------
}

Global_Part_Reload::~Global_Part_Reload()
{
  VEC_T::clean(epart); VEC_T::clean(npart);
}

// EOF
