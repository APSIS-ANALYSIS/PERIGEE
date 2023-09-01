#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
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

int main(int argc, char *argv[])
{
  SymmMatrix_3x3 smat; smat.gen_rand();
  Matrix_3x3 mat = smat.convert_to_full();

  Tensor4_3D ten = gen_T4_symm_id();
  SymmTensor4_3D sten = gen_ST4_symm_id();

  ten.gen_Ptilde( mat );
  sten.gen_Ptilde( smat );
  
  SymmMatrix_3x3 smat2; 
  std::vector<double> temp {}; 
  for(int ii=0; ii<100000; ++ii)
  {
    sten.gen_rand(-2.5, 12.2);
    for(int jj=0; jj<21; ++jj) temp.push_back( sten(jj) );
  }
  MATH_T::print_Histogram(temp);
  
  Matrix_3x3 mat2 = smat2.convert_to_full();

  //ten.add_SymmOutProduct(3.14159, mat, mat2);
  //sten.add_SymmOutProduct(3.14159, smat, smat2);

  /*
  for(int ii=0; ii<100; ++ii)
  {
    //Vector_3 vec1; vec1.gen_rand(); vec1.normalize();
    //Vector_3 vec2; vec2.gen_rand(); vec2.normalize();

    //ten.add_SymmOutProduct(1.1232789, vec1, vec2, vec1, vec2);
    //sten.add_SymmOutProduct(1.1232789, vec1, vec2);

    //ten = gen_Ptilde(mat2);
    //sten = gen_Ptilde(smat2);

    //sten.print_in_mat();

    sten.gen_rand(-1, 1);
    ten = sten.convert_to_full();

    SymmMatrix_3x3 smat3; smat3.gen_rand(-1, 1); 
    SymmMatrix_3x3 smat4; smat4.gen_rand(-1, 1); 

    Tensor4_3D PP; PP.gen_zero();
    PP.add_OutProduct(1.1523235904, smat3.convert_to_full(), smat4.convert_to_full());

    ten.TenPMult( PP );
    sten.TenPMult( PP );

    if( sten.is_identical(ten, 1.0e-15) ) std::cout<<"passed! \n";
    else std::cout<<"error. \n";
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }
  */  
  return EXIT_SUCCESS;
}

// EOF
