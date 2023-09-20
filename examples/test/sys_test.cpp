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
#include "DataVecStr.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"
#include "IIEN.hpp"
#include "Gmsh_FileIO.hpp"
#include "FEANode.hpp"

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::string part_file("part");
  int rank = 1;
  FEANode * fNode = new FEANode(part_file, rank);

  std::vector<int> id { 13, 12, 153, 82, 21, 19, 245 };

  std::vector<double> px(7), py(7), pz(7);

  fNode -> get_ctrlPts_xyz( 7, &id[0], &px[0], &py[0], &pz[0] );
  auto val = fNode -> get_ctrlPts_xyz( id );

  if( VEC_T::is_equal( val[0], px, 1.0e-27 ) ) std::cout<<"good x!\n";
  else std::cout<<"bad x \n";

  if( VEC_T::is_equal( val[1], py, 1.0e-27 ) ) std::cout<<"good y!\n";
  else std::cout<<"bad y \n";

  if( VEC_T::is_equal( val[2], pz, 1.0e-27 ) ) std::cout<<"good z!\n";
  else std::cout<<"bad z \n";

  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
