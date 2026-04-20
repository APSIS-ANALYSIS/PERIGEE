#include "Global_Part_Reload.hpp"
#include "HDF5_Reader.hpp"

Global_Part_Reload::Global_Part_Reload( int cpu_size, 
    int in_ncommon, bool isDualGraph,
    const std::string &element_part_name,
    const std::string &node_part_name )
{
  const std::string efName = element_part_name + ".h5";

  auto eh5r = SYS_T::make_unique<HDF5_Reader>(efName);
  
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

  // --------------------------------------------------------------------------
  const std::string nfName = node_part_name + ".h5";

  auto nh5r = SYS_T::make_unique<HDF5_Reader>(nfName);

  npart = nh5r -> read_intVector("/", "part");
  
  field_offset = nh5r -> read_intVector("/", "field_offset");

  SYS_T::print_fatal_if( cpu_size != cpusize, "Error: Global_Part_Reload cpu_size is incompatible with prior partition.\n" );
  SYS_T::print_fatal_if( in_ncommon != dual_edge_ncommon, "Error: Global_Part_Reload in_ncommon is incompatible with prior partition.\n" );
  SYS_T::print_fatal_if( isDualGraph != isDual, "Error: Global_Part_Reload isDualGraph is incompatible with prior partition.\n" );

  std::cout<<"=== Global partition loaded from "<<efName<<" and "<<nfName<<std::endl;
}

// EOF
