#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  SymmTensor2_3D cc; cc.gen_rand();
  
  Tensor2_3D F; F.gen_rand();

  Tensor2_3D out = F * cc * Ten2::transpose(F);

  cc.push_forward_stress(F);

  out.print_in_row();

  out += -cc.full();
  
  out.print_in_row();


  SymmTensor2_3D S_bar; S_bar.gen_rand();

  const auto S_iso_1 = STen2::gen_DEV_part(S_bar, cc);

  Tensor2_3D S_iso_2 = Ten4::gen_P(cc) * S_bar.full();

  S_iso_2.print_in_row();

  S_iso_2 += -1.0 * S_iso_1.full();

  S_iso_2.print_in_row();

  cc = STen2::gen_zero();
  cc.print_in_row();

  return EXIT_SUCCESS;
}

// EOF

