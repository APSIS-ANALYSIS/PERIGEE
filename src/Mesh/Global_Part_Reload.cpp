#include "Global_Part_Reload.hpp"

Global_Part_Reload::Global_Part_Reload( const int &cpu_size, 
    const int &in_ncommon, const bool &isDualGraph,
    const std::string &element_part_name,
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

  int cpusize = eh5r -> read_intScalar("/", "cpu_size");

  epart = eh5r -> read_intVector("/", "part");

  delete eh5r; H5Fclose( efile_id );

  // --------------------------------------------------------------------------
  const std::string nfName = node_part_name + ".h5";

  hid_t nfile_id = H5Fopen( nfName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * nh5r = new HDF5_Reader( nfile_id );

  npart = nh5r -> read_intVector("/", "part");

  delete nh5r; H5Fclose( nfile_id );
  // --------------------------------------------------------------------------

  SYS_T::print_fatal_if( cpu_size != cpusize, "Error: Global_Part_Reload cpu_size is incompatible with prior partition.\n" );
  SYS_T::print_fatal_if( in_ncommon != dual_edge_ncommon, "Error: Global_Part_Reload in_ncommon is incompatible with prior partition.\n" );
  SYS_T::print_fatal_if( isDualGraph != isDual, "Error: Global_Part_Reload isDualGraph is incompatible with prior partition.\n" );

  std::cout<<"=== Global partition loaded from "<<efName<<" and "<<nfName<<std::endl;
}

Global_Part_Reload::~Global_Part_Reload()
{
  VEC_T::clean(epart); VEC_T::clean(npart);
}

// EOF
