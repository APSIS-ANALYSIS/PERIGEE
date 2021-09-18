#include "HDF5_Reader.hpp"
#include "HDF5_Writer.hpp"

int main( int argc, char * argv[] )
{
  hid_t file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  int a = -1;
  h5w -> write_intScalar( "aa", a );

  delete h5w;
  H5Fclose(file_id);
  return EXIT_SUCCESS;
}

// EOF
