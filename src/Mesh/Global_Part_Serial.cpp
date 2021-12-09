#include "Global_Part_Serial.hpp"

Global_Part_Serial::Global_Part_Serial( const class IMesh * const &mesh,
   const std::string &element_part_name, const std::string &node_part_name )
: isMETIS(false), isDual(false), dual_edge_ncommon(0)
{
  const int nElem = mesh->get_nElem();
  const int nFunc = mesh->get_nFunc();
 
  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int ee=0; ee<nElem; ++ee) epart[ee] = 0;

  for(int nn=0; nn<nFunc; ++nn) npart[nn] = 0;

  const int cpu_size = 1;
  const bool isDualGraph = true;
  const int in_ncommon = 1;

  write_part_hdf5(element_part_name, epart, nElem, cpu_size, isDualGraph, in_ncommon, isMETIS );
  write_part_hdf5(node_part_name, npart, nFunc, cpu_size, isDualGraph, in_ncommon, isMETIS );

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_Serial::~Global_Part_Serial()
{
  delete [] epart; epart = nullptr;
  delete [] npart; npart = nullptr;
}


void Global_Part_Serial::write_part_hdf5( const std::string &fileName,
    const idx_t * const &part_in,
    const int &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName ); fName.append(".h5");

  // file creation
  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w -> write_intScalar("part_size", part_size);
  h5w -> write_intScalar("cpu_size", cpu_size);

  h5w->write_intScalar("part_isdual", ( part_isdual ? 1 : 0 ) );
  h5w->write_intScalar("in_ncommon", in_ncommon);

  h5w->write_intScalar("isMETIS", ( isMETIS ? 1 : 0 ) );
  h5w->write_intVector( "part", part_in, part_size );

  h5w->write_intScalar("isSerial", ( is_serial() ? 1 : 0 ) );

  delete h5w; H5Fclose(file_id);
}

// EOF
