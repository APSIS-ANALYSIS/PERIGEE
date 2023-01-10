#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "PLocAssem_Tet_VMS_NS_GenAlpha.hpp"
#include "PGAssem_NS_FEM.hpp"
#include "PTime_NS_Solver.hpp"

int main( int argc, char *argv[])
{
  std::string part_file("part");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);

  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  //pmat-> gen_id( pNode );

  pmat->gen_perm_bc(pNode, locnbc);

  delete pmat;
  PetscFinalize();
  return EXIT_SUCCESS;
}

//EOF
