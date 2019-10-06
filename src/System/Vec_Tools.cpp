#include "Vec_Tools.hpp"

void VEC_T::write_int_h5( const char * const &file_name, 
    const char * const &dataname,
    const std::vector<int> &value )
{
  std::string fName( file_name );
  fName.append(".h5");

  hid_t file_id;
  file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer( file_id );

  h5writer->write_intVector(dataname, value);

  delete h5writer;
  H5Fclose(file_id);
}

void VEC_T::read_int_h5( const char * const &file_name,
    const char * const &groupname,
    const char * const &dataname,
    std::vector<int> &value )
{
  std::string fName( file_name );
  fName.append(".h5");

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5reader = new HDF5_Reader( file_id );

  hid_t drank;
  hsize_t * ddims;
  int * intarray;

  h5reader->read_intArray( groupname, dataname, drank, ddims, intarray );
  
  if(drank > 1)
  {
    std::stringstream ss; ss<<"Error: "<<dataname<<" at "<<groupname<<" "
      <<file_name<<".h5 is not an one-dimensional array. \n";
    SYS_T::commPrint(ss.str().c_str()); 
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  const int length = ddims[0];

  VEC_T::fillArray(value, intarray, length);

  delete [] ddims; delete [] intarray;
  delete h5reader;
  H5Fclose( file_id );
}

// EOF
