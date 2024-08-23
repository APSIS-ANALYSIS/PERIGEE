#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "ElemBC_3D_turbulence_wall_model.hpp"
#include "EBC_Partition_outflow.hpp"
#include "EBC_Partition_turbulence_wall_model.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include <iomanip>
#include "PETSc_Tools.hpp"

void test(const double * const &vv, const unsigned int len)
{
  for(unsigned int ii=0; ii<len; ++ii)
    std::cout<<vv[ii]<<'\t';
  std::cout<<"End of vector \n";
}

int main(int argc, char *argv[])
{
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  Vec GG;
  VecCreate(PETSC_COMM_WORLD, &GG);
  VecSetFromOptions(GG);
  VecSetSizes(GG, 10, PETSC_DECIDE);

  PetscInt idx[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  PetscScalar val[] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1};

  VecSetValues(GG, 10, idx, val, INSERT_VALUES);

  VecAssemblyBegin(GG);
  VecAssemblyEnd(GG);

  VecView(GG, PETSC_VIEWER_STDOUT_WORLD);

  PetscInt * idx_from = new PetscInt[4];
  idx_from[0] = 1;
  idx_from[1] = -1;
  idx_from[2] = 5;
  idx_from[3] = 7;
  VecSetOption(GG, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  double * sca = new double [4];
  PETSc_T::Scatter(GG, idx_from, 4, sca);

  for(int ii=0; ii<4; ++ii)
    std::cout << idx_from[ii] << ':' << sca[ii] << std::endl;

  VecDestroy(&GG);

  delete [] sca;

  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
