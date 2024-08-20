#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  Vector_3 a; a.gen_rand();

  Tensor2_3D B = Ten2::gen_zero();

  B.gen_outprod(a);

  auto C = STen2::gen_dyad(a);

  B.print_in_row();
  
  B -= C.full();

  B.print_in_row();

  return EXIT_SUCCESS;
}

// EOF

