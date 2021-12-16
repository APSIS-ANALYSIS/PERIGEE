#include "Global_Part_Serial.hpp"

Global_Part_Serial::Global_Part_Serial( const class IMesh * const &mesh,
   const std::string &element_part_name, const std::string &node_part_name )
: isMETIS(false), isDual(false), dual_edge_ncommon(0)
{
  // This is a partition for a single mesh (field)
  field_offset.resize(1);
  field_offset[0] = 0;

  const int nElem = mesh->get_nElem();
  const int nFunc = mesh->get_nFunc();
 
  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int ee=0; ee<nElem; ++ee) epart[ee] = 0;

  for(int nn=0; nn<nFunc; ++nn) npart[nn] = 0;

  const int cpu_size = 1;

  write_part_hdf5(element_part_name, epart, nElem, cpu_size );
  write_part_hdf5(node_part_name, npart, nFunc, cpu_size );

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_Serial::Global_Part_Serial( const int &num_fields,
    const std::vector<IMesh const *> &mesh_list,
    const std::string &element_part_name, const std::string &node_part_name )
: isMETIS( false ), isDual( false ), dual_edge_ncommon( 0 )
{
  if(num_fields != static_cast<int>( mesh_list.size() ) )
  {
    std::cerr<<"ERROR: input num_fields is incompatible with mesh list.\n";
    exit(1);
  }

  field_offset.resize( num_fields );
  field_offset[0] = 0;

  for(int ii=1; ii<num_fields; ++ii)
    field_offset[ii] = field_offset[ii-1] + mesh_list[ii] -> get_nFunc();

  const idx_t nElem = mesh_list[0]->get_nElem();
  idx_t nFunc = 0;
  idx_t nLocBas = 0;

  for(int ii=0; ii<num_fields; ++ii)
  {
    if( nElem != static_cast<idx_t>( mesh_list[ii]->get_nElem() ) )
    {
      std::cerr<<"ERROR: mesh list objects are incompatible with nElem.\n";
      exit(1);
    }
    nFunc   += mesh_list[ii] -> get_nFunc();
    nLocBas += mesh_list[ii] -> get_nLocBas();
  }

  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int ee=0; ee<nElem; ++ee) epart[ee] = 0;

  for(int nn=0; nn<nFunc; ++nn) npart[nn] = 0;

  int cpu_size = 1;

  write_part_hdf5(element_part_name, epart, nElem, cpu_size );
  write_part_hdf5(node_part_name, npart, nFunc, cpu_size );

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_Serial::~Global_Part_Serial()
{
  delete [] epart; epart = nullptr;
  delete [] npart; npart = nullptr;
}

void Global_Part_Serial::write_part_hdf5( const std::string &fileName,
    const idx_t * const &part_in,
    const int &part_size, const int &cpu_size ) const
{
  std::string fName( fileName ); fName.append(".h5");

  // file creation
  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w -> write_intScalar("part_size", part_size);
  h5w -> write_intScalar("cpu_size", cpu_size);

  h5w->write_intScalar("part_isdual", ( isDual ? 1 : 0 ) );
  h5w->write_intScalar("in_ncommon", dual_edge_ncommon);

  h5w->write_intScalar("isMETIS", ( isMETIS ? 1 : 0 ) );
  h5w->write_intVector("part", part_in, part_size );

  h5w->write_intVector("field_offset", field_offset );

  h5w->write_intScalar("isSerial", ( is_serial() ? 1 : 0 ) );

  delete h5w; H5Fclose(file_id);
}

// EOF
