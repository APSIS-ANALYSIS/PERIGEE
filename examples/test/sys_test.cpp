#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor2_3D.hpp"
#include "SymmTensor2_3D.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "VTK_Tools.hpp"
#include "NodalBC.hpp"
#include "Gmsh_FileIO.hpp"

int main(int argc, char *argv[])
{
  double mean = 1.5, std = 5.23;

  std::vector<double> val1, val2;
  MATH_T::gen_Gaussian(10000, mean, std, val1);

  val2.clear();
  for(int ii=0; ii<10000; ++ii)
    val2.push_back( MATH_T::gen_double_rand_normal(mean, std) );

  MATH_T::print_Histogram(val1);
  std::cout<<std::endl;
  MATH_T::print_Histogram(val2);

  return EXIT_SUCCESS;
}

// EOF
