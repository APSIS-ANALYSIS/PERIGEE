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

int main(int argc, char *argv[])
{

  HEX_T::Hex8 * test_hex = new HEX_T::Hex8();  
  test_hex -> reset(8, 10, 1, 6, 21, 22, 5, 7);

  std::cout << "face_id:" << test_hex->get_face_id(10, 1, 22, 5) << std::endl;
  std::cout << "face_id:" << test_hex->get_face_id(22, 10, 1, 5) << std::endl;
  //std::cout << "face_id:" << test_hex->get_face_id(7, 10, 1, 5) << std::endl;
  std::cout << "face_id:" << test_hex->get_face_id(23, 10, 1, 5) << std::endl;

  return EXIT_SUCCESS;
}

// EOF
