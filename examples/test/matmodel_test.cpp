#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  Tensor4_3D AA = Ten4::gen_zero();
  auto BB = AA;

  SymmTensor2_3D cc; cc.gen_rand();
  Tensor2_3D dd = cc.full();
  
  SymmTensor2_3D ee; ee.gen_rand();
  Tensor2_3D ff = ee.full();

  AA.add_OutProduct(3.115, dd, ff );
  BB.add_OutProduct(3.115, cc, ee );

  AA.print();
  AA.AXPY(-1.0, BB);

  AA.print();

  return EXIT_SUCCESS;
}

// EOF

