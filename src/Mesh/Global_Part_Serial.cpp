#include "Global_Part_Serial.hpp"

Global_Part_Serial::Global_Part_Serial( const int &in_nelem, const int &in_nfunc,
   const std::string &element_part_name, const std::string &node_part_name )
{
  // This is a partition for a single mesh (field)
  field_offset.resize(1);
  field_offset[0] = 0;

  const int nElem = in_nelem;
  const int nFunc = in_nfunc;
 
  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int ee=0; ee<nElem; ++ee) epart[ee] = 0;

  for(int nn=0; nn<nFunc; ++nn) npart[nn] = 0;

  write_part_hdf5(element_part_name, epart, nElem );
  write_part_hdf5(node_part_name, npart, nFunc );

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_Serial::Global_Part_Serial( const int &num_fields,
    const std::vector<int> &nelem_list,
    const std::vector<int> &nfunc_list,
    const std::string &element_part_name, const std::string &node_part_name )
{
  if(num_fields != static_cast<int>( nelem_list.size() ) )
  {
    std::cerr<<"ERROR: input num_fields is incompatible with nelem list.\n";
    exit(1);
  }

  if(num_fields != static_cast<int>( nfunc_list.size() ) )
  {
    std::cerr<<"ERROR: input num_fields is incompatible with nfunc list.\n";
    exit(1);
  }

  field_offset.resize( num_fields );
  field_offset[0] = 0;

  for(int ii=1; ii<num_fields; ++ii)
    field_offset[ii] = field_offset[ii-1] + nfunc_list[ii-1];

  const idx_t nElem = nelem_list[0];
  idx_t nFunc = 0;

  for(int ii=0; ii<num_fields; ++ii)
  {
    if( nElem != static_cast<idx_t>( nelem_list[ii] ) )
    {
      std::cerr<<"ERROR: nelem list objects are incompatible with nElem.\n";
      exit(1);
    }
    nFunc += nfunc_list[ii];
  }

  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int ee=0; ee<nElem; ++ee) epart[ee] = 0;

  for(int nn=0; nn<nFunc; ++nn) npart[nn] = 0;

  write_part_hdf5(element_part_name, epart, nElem);
  write_part_hdf5(node_part_name, npart, nFunc);

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_Serial::~Global_Part_Serial()
{
  delete [] epart; epart = nullptr; delete [] npart; npart = nullptr;
}

void Global_Part_Serial::write_part_hdf5( const std::string &fileName,
    const idx_t * const &part_in, const int &part_size ) const
{
  const std::string fName = fileName + ".h5";

  // file creation
  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  auto h5w = SYS_T::make_unique<HDF5_Writer>(file_id);

  h5w->write_intScalar("part_size", part_size);
  h5w->write_intScalar("cpu_size", 1);
  h5w->write_intVector("part", part_in, part_size );
  h5w->write_intVector("field_offset", field_offset );

  h5w->write_intScalar("isMETIS", ( get_isMETIS() ? 1 : 0 ) );
  h5w->write_intScalar("part_isdual", ( get_isDual() ? 1 : 0 ) );
  h5w->write_intScalar("in_ncommon", get_dual_edge_ncommon() );
  
  h5w->write_intScalar("isSerial", ( is_serial() ? 1 : 0 ) );

  H5Fclose(file_id);
}

// EOF
