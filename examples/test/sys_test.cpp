#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "HDF5_Reader.hpp"

int main( int argc, char * argv[] )
{
  std::vector<std::string> a;
  a.push_back("hello a\n");
  a.push_back("cmame \t");
  a.push_back("july 31th \n");
  a.push_back("aug 31th \n");

  for(auto &out : a) out = "hello\n";

  //for(auto out : a) std::cout<<out;

  std::string fName = "part_p00000.h5";

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  if( h5r -> check_data("/ebc_wall/ebcid_0/thicks") )
    std::cout<<"prestress exist\n";
  else
    std::cout<<"prestress does not exist\n";

  H5Fclose( file_id ); delete h5r;

  return EXIT_SUCCESS;
}

// EOF
