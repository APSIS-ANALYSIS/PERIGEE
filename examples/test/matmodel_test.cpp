#include "Sys_Tools.hpp"
#include "Vector_3.hpp"

int main( int argc, char * argv[] )
{
  const std::array<double, 3> aa {{ 1.2, 3.1, -1.2 }};

  const Vector_3 bb(aa);

  bb.print();

  return EXIT_SUCCESS;
}

// EOF

