#include "HDF5_Reader.hpp"
#include "HDF5_Writer.hpp"

int main( int argc, char * argv[] )
{
  hid_t file_id = H5Fcreate( "test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  int a = -1;
  h5w -> write_intScalar( "aa", a );

  /*            test write_Vector_3 and read_Vector_3 functions           */
  Vector_3 vv;

  vv.gen_rand();

  h5w -> write_Vector_3( "vv", vv );

  hid_t group_id = H5Gcreate( file_id, "/GROUP_TEST",
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  
  h5w -> write_Vector_3( group_id, "vv", vv );

  H5Gclose( group_id );

  delete h5w;
  H5Fclose( file_id );

  hid_t file_id_2 = H5Fopen( "test.h5", H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id_2 );

  Vector_3 rr = h5r -> read_Vector_3("/GROUP_TEST", "vv");

  delete h5r; H5Fclose(file_id_2);

  rr.print();

  rr -= vv;

  rr.print();

  /*            test write_Matrix_3x3 and read_Matrix_3x3 functions           */

  hid_t file2 = H5Fcreate( "test2.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  
  HDF5_Writer * h5w2 = new HDF5_Writer( file2 );
  
  Matrix_3x3 mm;
  
  mm.gen_rand();
  
  h5w2 -> write_Matrix_3x3( "mm", mm);

  hid_t group2 = H5Gcreate( file2, "/GROUP2_TEST",
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  h5w2 -> write_Matrix_3x3( group2, "mm", mm );

  H5Gclose( group2 );

  delete h5w2;
  H5Fclose( file2 );

  hid_t file3 = H5Fopen( "test2.h5", H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r2 = new HDF5_Reader( file3 );

  Matrix_3x3 rmm = h5r2 -> read_Matrix_3x3( "/GROUP2_TEST", "mm" );

  delete h5r2; H5Fclose( file3 );

  rmm.print();

  rmm -= mm;

  rmm.print();

  return EXIT_SUCCESS;
}

// EOF
