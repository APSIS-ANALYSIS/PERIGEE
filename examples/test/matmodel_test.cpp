#include "ANL_Tools.hpp"

int main( int argc, char * argv[] )
{
  std::string part_file("part");

  int size = ANL_T::get_cpu_size(part_file, 0);

  for(int ii=0; ii<size; ++ii)
  {
    int rank = ANL_T::get_cpu_rank(part_file, ii);
    size = ANL_T::get_cpu_size(part_file, ii);
    std::cout<<rank<<'\t'<<size<<'\n';
  }

  return 0;
}

// EOF

