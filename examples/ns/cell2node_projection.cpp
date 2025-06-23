#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_NS.hpp"
#include "APart_Node.hpp"
#include "FEANode.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_NBC.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"
#include "ANL_Tools.hpp"
#include "VisDataPrep_NS.hpp"
#include "VTK_Writer_NS.hpp"

// Zero_Tangent_Rezidual
void ZTR(PetscScalar * T, PetscScalar * R, const int &vecsize);

int main(int argc, char *argv[])
{
  using namespace std;
  int nqpv = 5;
  int dof_mat = 4;
  std::string part_file("part");
  std::string sol_bName("SOL_"); // base name of the solution file

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  // PetscInitialize(&argc, &argv, NULL, NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  // MPI_Init(&argc, &argv);
  // int rank, size;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout<< "===> Rank = " << rank << endl;
  cout<<"===> Size = " << size << endl;

  SYS_T::commPrint("===> Setup analysis tools.\n");
  // Local sub-domain's nodal indices
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  // Local sub-domain's IEN array
  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);

  // Local sub-domain's nodal bc
  auto locnbc = SYS_T::make_unique<ALocal_NBC>(part_file, rank);

  SYS_T::commPrint("===> Setup element and quadrature rule.\n");
  // Element and Quadrature rules
  auto elemType = ANL_T::get_elemType(part_file, rank);
  auto elementv = ElementFactory::createVolElement(elemType, nqpv);
  auto quadv = QuadPtsFactory::createVolQuadrature(elemType, nqpv);
  const int nLocBas = elementv->get_nLocBas();
  const int nqp = quadv->get_num_quadPts();

  // Global matrix and vectors
  SYS_T::commPrint("===> Initialize global matrix and vectors.\n");
  Mat L2_proj_mat;
  Vec L2_proj_rhs;

  const int nlocrow = dof_mat * pNode->get_nlocalnode();
  
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, 
    PETSC_DETERMINE, PETSC_DETERMINE,
    PETSC_DECIDE, NULL, PETSC_DECIDE, NULL, &L2_proj_mat);

  MatZeroEntries(L2_proj_mat);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_rhs);
  VecSetSizes(L2_proj_rhs, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(L2_proj_rhs);
  VecSet(L2_proj_rhs, 0.0);
  VecSetOption(L2_proj_rhs, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  std::unique_ptr<PDNSolution> sol =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), 0 );
  
  // Local matrix and vectors
  SYS_T::commPrint("===> Initialize local matrix and vector.\n");
  const int vec_size = nLocBas * dof_mat;
  auto Tangent = new PetscScalar[vec_size * vec_size];
  auto Residual = new PetscScalar[vec_size];

  ZTR(Tangent, Residual, vec_size);
  
  // Solver
  SYS_T::commPrint("===> Initialize solver.\n");
  auto lsolver_L2_proj = SYS_T::make_unique<PLinear_Solver_PETSc>(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

  KSPSetType(lsolver_L2_proj->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(lsolver_L2_proj->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(lsolver_L2_proj->ksp, 500);

  PC preproc; lsolver_L2_proj->GetPC(&preproc);
  PCSetType( preproc, PCHYPRE );
  PCHYPRESetType( preproc, "boomeramg" );
  
  // Loop and assemble
  SYS_T::commPrint("===> Assembling......\n");
  auto ctrl_x = new double [nLocBas];
  auto ctrl_y = new double [nLocBas];
  auto ctrl_z = new double [nLocBas];
  auto IEN_e = new int [nLocBas];
  auto row_index = new PetscInt [vec_size];

  std::ostringstream time_index;

  for(int tt = 500; tt < 18000; tt += 500)
  {
    SYS_T::commPrint("  Time = %d\n", tt);
    MatZeroEntries(L2_proj_mat);
    VecSet(L2_proj_rhs, 0.0);

    // Control points' xyz coordinates
    auto fNode = SYS_T::make_unique<FEANode>(part_file, rank, tt);

    // Local sub-domain's element indices
    auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank, tt);

    const int nelem = locElem->get_nlocalele();

    for(int ee=0; ee<nelem; ++ee)
    {
      // Cell data
      const double pp = locElem->get_local_p(ee);
      const double uu = locElem->get_local_u(ee);
      const double vv = locElem->get_local_v(ee);
      const double ww = locElem->get_local_w(ee);

      locIEN->get_LIEN(ee, IEN_e);

      fNode->get_ctrlPts_xyz(nLocBas, IEN_e, ctrl_x, ctrl_y, ctrl_z);

      elementv->buildBasis(quadv.get(), ctrl_x, ctrl_y, ctrl_z);

      ZTR(Tangent, Residual, vec_size);

      for(int qua=0; qua<nqp; ++qua)
      {
        std::vector<double> R = elementv->get_R(qua);

        const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

        for(int A=0; A<nLocBas; ++A)
        {
          const double NA = R[A];
          
          const int A4 = 4 * A;

          Residual[A4] += gwts * NA * pp;
          Residual[A4+1] += gwts * NA * uu;
          Residual[A4+2] += gwts * NA * vv;
          Residual[A4+3] += gwts * NA * ww;

          for(int B=0; B<nLocBas; ++B)
          {
            const double NB = R[B];

            Tangent[16*nLocBas*A+4*B] += gwts * NA * NB;

            Tangent[4*nLocBas*(4*A+1)+4*B+1] += gwts * NA * NB;

            Tangent[4*nLocBas*(4*A+2)+4*B+2] += gwts * NA * NB;

            Tangent[4*nLocBas*(4*A+3)+4*B+3] += gwts * NA * NB;
          }
        }
      }

      for(int ii=0; ii<nLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          row_index[dof_mat*ii+mm] = dof_mat * locnbc -> get_LID(mm, IEN_e[ii]) + mm;
      }

      MatSetValues(L2_proj_mat, vec_size, row_index, vec_size, row_index,
        Tangent, ADD_VALUES);

      VecSetValues(L2_proj_rhs, vec_size, row_index,
        Residual, ADD_VALUES);
    }

    MatAssemblyBegin(L2_proj_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(L2_proj_mat, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(L2_proj_rhs);
    VecAssemblyEnd(L2_proj_rhs);

    // Solve
    SYS_T::commPrint("===> Solving......\n");
    lsolver_L2_proj->Solve(L2_proj_mat, L2_proj_rhs, sol.get(), true);
    SYS_T::commPrint("\n");

    std::string name_to_write(sol_bName);

    time_index.str("");
    time_index << 900000000 + tt;

    name_to_write.append(time_index.str());

    // Write solution
    sol->WriteBinary(name_to_write);
  }

  // Finalization
  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z;
  delete [] IEN_e; delete [] Tangent; delete [] Residual;

  VecDestroy(&L2_proj_rhs);
  MatDestroy(&L2_proj_mat);

  SYS_T::commPrint("===> Finished!\n");

  PetscFinalize();

  // MPI_Finalize();

  return EXIT_SUCCESS;
}

void ZTR(PetscScalar * T, PetscScalar * R, const int &vecsize)
{
  for(int ii=0; ii<vecsize; ++ii)
    R[ii] = 0.0;

  for(int ii=0; ii<vecsize * vecsize; ++ii)
    T[ii] = 0.0;
}