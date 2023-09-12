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
#include "DataVecStr.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"
#include "IIEN.hpp"
#include "Gmsh_FileIO.hpp"

int main(int argc, char *argv[])
{
  std::string gmshFile = "two_cube.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  // GIO -> print_info();
  GIO -> update_quadratic_hex_IEN(0);
  GIO -> update_quadratic_hex_IEN(1);

  GIO -> write_quadratic_sur_vtu("ftop",0,0,true);
  GIO -> write_quadratic_sur_vtu("fbot",1,0,true);
  GIO -> write_quadratic_sur_vtu("fwall",2,0,true);

  GIO -> write_quadratic_sur_vtu("fswall",2,1,true);
  GIO -> write_quadratic_sur_vtu("sbot",3,1);
  GIO -> write_quadratic_sur_vtu("stop",4,1);
  GIO -> write_quadratic_sur_vtu("swall",5,1);

  GIO -> write_each_vtu();

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
