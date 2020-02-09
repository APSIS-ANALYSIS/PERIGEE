#include "Global_Part_Serial.hpp"

Global_Part_Serial::Global_Part_Serial( const class IMesh * const &mesh,
   const char * const &element_part_name, const char * const &node_part_name )
{
  const int nElem = mesh->get_nElem();
  const int nFunc = mesh->get_nFunc();
 
  isMETIS = false;
  isDual = false;
  dual_edge_ncommon = 0;

  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];

  for(int e = 0; e<nElem; ++e) epart[e] = 0;

  for(int n=0; n<nFunc; ++n) npart[n] = 0;

  int cpu_size = 1;
  bool isDualGraph = true;
  int in_ncommon = 1;

  write_part_hdf5(element_part_name, epart, nElem, cpu_size, isDualGraph, in_ncommon, isMETIS );
  write_part_hdf5(node_part_name, npart, nFunc, cpu_size, isDualGraph, in_ncommon, isMETIS );

  std::cout<<"\n=== Global partition generated. \n";
}

Global_Part_Serial::~Global_Part_Serial()
{
  delete [] epart; epart = NULL;
  delete [] npart; npart = NULL;
  std::cout<<"-- Global_Part_Serial deleted. \n";
}


inline idx_t Global_Part_Serial::get_epart( int e ) const
{
  return epart[e];
}


inline idx_t Global_Part_Serial::get_npart( int n ) const
{
  return npart[n];
}


void Global_Part_Serial::write_part_hdf5( const char * const &fileName,
    const idx_t * const &part_in,
    const int &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName );
  fName.append(".h5");

  hid_t file_id, dataspace_id_n, dataspace_id_1;
  hsize_t dim_n[1], dim_1[1];

  // file creation
  file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // dataspace creation
  dim_n[0] = part_size; dim_1[0] = 1;

  dataspace_id_n = H5Screate_simple(1, dim_n, NULL);
  dataspace_id_1 = H5Screate_simple(1, dim_1, NULL);

  // dataset
  hid_t setid_part_size, setid_cpu_size, setid_part_isdual, setid_in_ncommon;
  hid_t setid_isMETIS, setid_part;

  setid_part_size = H5Dcreate( file_id, "part_size", H5T_NATIVE_INT, dataspace_id_1,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_cpu_size = H5Dcreate( file_id, "cpu_size", H5T_NATIVE_INT, dataspace_id_1,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_part_isdual = H5Dcreate( file_id, "part_isdual", H5T_NATIVE_INT, dataspace_id_1,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_in_ncommon = H5Dcreate( file_id, "in_ncommon", H5T_NATIVE_INT, dataspace_id_1,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_isMETIS = H5Dcreate( file_id, "isMETIS", H5T_NATIVE_INT, dataspace_id_1,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_part = H5Dcreate( file_id, "part", H5T_NATIVE_INT, dataspace_id_n,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


  H5Dwrite( setid_part_size, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &part_size);
  H5Dwrite( setid_cpu_size, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &cpu_size);
  
  int intbool;
  if(part_isdual)
   intbool = 1;
  else
   intbool = 0; 
  H5Dwrite( setid_part_isdual, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &intbool);

  H5Dwrite( setid_in_ncommon, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &in_ncommon);

  if(isMETIS)
   intbool = 1;
  else
   intbool = 0; 
  H5Dwrite( setid_isMETIS, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &intbool);

  H5Dwrite( setid_part, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
     part_in );

  H5Dclose( setid_part_size );   H5Dclose( setid_cpu_size );
  H5Dclose( setid_part_isdual ); H5Dclose( setid_in_ncommon );
  H5Dclose( setid_isMETIS );     H5Dclose( setid_part );
  H5Sclose( dataspace_id_n );
  H5Sclose( dataspace_id_1 );
  H5Fclose(file_id);
}

// EOF
