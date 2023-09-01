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
  
  SymmMatrix_3x3 smat2; smat2.gen_rand(); 
  Matrix_3x3 mat2 = smat2.convert_to_full();

  //ten.add_SymmOutProduct(3.14159, mat, mat2);
  //sten.add_SymmOutProduct(3.14159, smat, smat2);

  for(int ii=0; ii<100; ++ii)
  {
    sten.gen_rand();
    ten = sten.convert_to_full();
    
    const double rval = MATH_T::gen_randomD_closed(-1.11, 1.23);

    /*
    Vector_3 vec1; vec1.gen_rand(); vec1.normalize();
    Vector_3 vec2; vec2.gen_rand(); vec2.normalize();


    ten.add_SymmOutProduct(rval, vec1, vec2, vec1, vec2);
    sten.add_SymmOutProduct(rval, vec1, vec2);
    */

    SymmMatrix_3x3 smat2; smat2.gen_rand();
    Matrix_3x3 mat2 = smat2.convert_to_full();

    ten.add_SymmOutProduct(rval, mat, mat2);
    sten.add_SymmOutProduct(rval, smat, smat2);

    //ten = gen_Ptilde(mat2);
    //sten = gen_Ptilde(smat2);

    //sten.print_in_mat();

    /*
       sten.gen_rand(-1, 1);
       ten = sten.convert_to_full();

       SymmMatrix_3x3 smat3; smat3.gen_rand(-1, 1); 
       SymmMatrix_3x3 smat4; smat4.gen_rand(-1, 1); 

       Tensor4_3D PP; PP.gen_zero();
       PP.add_OutProduct(1.1523235904, smat3.convert_to_full(), smat4.convert_to_full());

       ten.TenPMult( PP );
       sten.TenPMult( PP );
       */

    if( sten.is_identical(ten, 2.0e-15) ) std::cout<<"passed! \n";
    else std::cout<<"error. \n";
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  SymmTensor4_3D aaa = gen_ST4_symm_id();
  
  SymmTensor4_3D bbb; bbb.gen_symm_id();
  
  aaa -= bbb;

  aaa.print_in_mat();

  Tensor4_3D ta = gen_T4_symm_id();
  Tensor4_3D tb; tb.gen_symm_id();

  if( ta.is_identical(tb, 1.0e-17) ) std::cout<<"symm_id is good! \n";
  else std::cout<<"symm_id is bad. \n";

  ta -= tb;

  ta.print_in_mat();

  return EXIT_SUCCESS;
}

// EOF
