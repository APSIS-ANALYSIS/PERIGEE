#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  SymmTensor2_3D a;
  a.gen_rand();

  auto v = a.to_std_vector();

  VEC_T::print(v);

  return EXIT_SUCCESS;
}

// EOF

