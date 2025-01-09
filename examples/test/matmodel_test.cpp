#include "APart_Basic_Info.hpp"

int main( int argc, char * argv[] )
{
  std::string part_file("part");

  int size = APart_Basic_Info::get_cpu_size(part_file, 0);

  for(int ii=0; ii<size; ++ii)
  {
    int rank = APart_Basic_Info::get_cpu_rank(part_file, ii);
    size = APart_Basic_Info::get_cpu_size(part_file, ii);
    std::cout<<rank<<'\t'<<size<<'\n';
  }

  return 0;
}

// EOF

