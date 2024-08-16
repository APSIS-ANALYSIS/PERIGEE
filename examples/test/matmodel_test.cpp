#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  SymmTensor2_3D cc; cc.gen_rand();
  
  Tensor2_3D F; F.gen_rand(-10, 10);

  Tensor2_3D out = F * cc * Ten2::transpose(F);

  cc.push_forward_stress(F);

  out.print_in_row();

  out += -cc.full();
  
  out.print_in_row();

  return EXIT_SUCCESS;
}

// EOF

