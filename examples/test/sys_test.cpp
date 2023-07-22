#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"

int main(int argc, char *argv[])
{
  const double tol = 1e-13;
  const int max_itr = 1000;
  for(int i = 0; i < max_itr; ++i)
  {
    SymmMatrix_3x3 symm_mat_1;
    symm_mat_1.gen_rand();

    SymmMatrix_3x3 symm_mat_2;
    symm_mat_2.gen_rand();

    Matrix_3x3 mat_2 = symm_mat_2.convert_to_full();

    double symm_contraction = symm_mat_1.MatContraction( symm_mat_2 );
    double contraction = symm_mat_1.MatContraction( mat_2 );
    
    std::cout<<std::setprecision(8)<<"step "<<i<<" : "<<symm_contraction<<" "<<contraction<<std::endl;
    if (std::abs(symm_contraction - contraction) > tol)
    {
      std::cout<<"Warning: unsafe with tol: "<<tol<<std::endl;
      break;
    }
    usleep(1000000);
  }

  return EXIT_SUCCESS;
}

// EOF
