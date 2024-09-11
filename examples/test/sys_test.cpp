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
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_vol_ST91.hpp"
#include "MaterialModel_vol_M94.hpp"
#include "MaterialModel_ich_GOH06.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "MaterialModel_ich_StVenant_Kirchhoff.hpp"
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
  VecSetSizes(GG, 100, PETSC_DECIDE);

  PetscInt * idx = new PetscInt[100];
  PetscScalar * val = new PetscScalar[100];

  const int rank = SYS_T::get_MPI_rank();
  for(int ii=0; ii<100; ++ii)
  {
    idx[ii] = 100 * rank + ii;
    val[ii] = double(100 * rank + ii) + 0.1;
  }

  VecSetValues(GG, 100, idx, val, INSERT_VALUES);

  VecAssemblyBegin(GG);
  VecAssemblyEnd(GG);

  PetscInt * idx_from = new PetscInt[4];
  idx_from[0] = 1;
  idx_from[1] = 11;
  idx_from[2] = 21;
  idx_from[3] = 31;
  double * sca = new double [4];


  for(int iter=0; iter< 5-rank; ++iter)
  {
    PETSc_T::Scatter(GG, idx_from, 4, sca);

    if(rank == 0)
      for(int ii=0; ii<4; ++ii)
        std::cout << idx_from[ii] << ':' << sca[ii] << std::endl;
  }

  // Synchronize
  for(int iter=0; iter<rank; ++iter)
  {
    PETSc_T::Scatter(GG, idx_from, 4, sca);

    if(rank == 0)
      for(int ii=0; ii<4; ++ii)
        std::cout << idx_from[ii] << ':' << sca[ii] << std::endl;
  }

  VecDestroy(&GG);

  delete [] idx;
  delete [] idx_from;
  delete [] val;
  delete [] sca;

  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
