#include "PLinear_Solver_2x2Block.hpp"

PLinear_Solver_2x2Block::PLinear_Solver_2x2Block( 
    const double &rtol0, const double &atol0, 
    const double &dtol0, const int &maxit0,
    const double &rtol1, const double &atol1, 
    const double &dtol1, const int &maxit1,
    const IPGAssem_2x2Block * const &gAssem_ptr )
{
  solver_0 = new PLinear_Solver_PETSc( rtol0, atol0, dtol0, maxit0, 
      "ls0_", "pc0_" );

  solver_1 = new PLinear_Solver_PETSc( rtol1, atol1, dtol1, maxit1,
      "ls1_", "pc1_" );

  // Allocate G0 and G1
  VecDuplicate( gAssem_ptr->G_0, &G0 );
  VecDuplicate( gAssem_ptr->G_1, &G1 );

  // Allocate D0 and D1
  VecDuplicate( gAssem_ptr->G_0, &D0 );
  VecDuplicate( gAssem_ptr->G_1, &D1 );
}


PLinear_Solver_2x2Block::~PLinear_Solver_2x2Block()
{
  delete solver_0; delete solver_1;
  VecDestroy(&G0); VecDestroy(&G1);
  VecDestroy(&D0); VecDestroy(&D1);
}


void PLinear_Solver_2x2Block::SymmJacobi_Scale( 
    IPGAssem_2x2Block * const &gAssem_ptr )
{
  // Obtain diagonal entries
  PETSc_T::GetDiagonal(gAssem_ptr->K_00, D0);
  PETSc_T::GetDiagonal(gAssem_ptr->K_11, D1);

  // Power to -1/2
  PETSc_T::MinusSqrtVec( D0 );
  PETSc_T::MinusSqrtVec( D1 );

  // Upage K matrices
  PETSc_T::DiagonalScale(gAssem_ptr->K_00, D0, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_01, D0, D1);
  PETSc_T::DiagonalScale(gAssem_ptr->K_10, D1, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_11, D1, D1);

  // Update G vectors
  PETSc_T::DiagonalScale( D0, gAssem_ptr->G_0 );
  PETSc_T::DiagonalScale( D1, gAssem_ptr->G_1 );
}


void PLinear_Solver_2x2Block::SymmJacobi_Restore(
    IPGAssem_2x2Block * const &gAssem_ptr )
{
  // Power to -1
  PETSc_T::InvAbsVec( D0 );
  PETSc_T::InvAbsVec( D1 );

  // Upage K matrices
  PETSc_T::DiagonalScale(gAssem_ptr->K_00, D0, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_01, D0, D1);
  PETSc_T::DiagonalScale(gAssem_ptr->K_10, D1, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_11, D1, D1);

  // Update G vectors
  PETSc_T::DiagonalScale( D0, gAssem_ptr->G_0 );
  PETSc_T::DiagonalScale( D1, gAssem_ptr->G_1 );
}


void PLinear_Solver_2x2Block::Solve_ExaFac( 
    const IPGAssem_2x2Block * const &gAssem_ptr,
    PDNSolution * const &sol_0, PDNSolution * const &sol_1 )
{
  // Create the Schur complement S
  MatCreateSchurComplement(gAssem_ptr->K_11, gAssem_ptr->K_11,
      gAssem_ptr->K_10, gAssem_ptr->K_01, gAssem_ptr->K_00, &S);

  // Assign operators
  solver_0 -> SetOperator( S );
  solver_1 -> SetOperator( gAssem_ptr->K_11 );

  // Setup external solvers
  KSPSetType(solver_0 -> ksp, KSPCG);
  KSPSetType(solver_1 -> ksp, KSPCG);

  PC pc0; solver_0 -> GetPC(&pc0); PCSetType(pc0, PCNONE);
  PC pc1; solver_1 -> GetPC(&pc1); PCSetType(pc1, PCJACOBI);

  // Setup internal solver
  KSP ksp_s;
  MatSchurComplementGetKSP(S, &ksp_s);
  const double ksp_s_rtol = 1.0e-2;
  const double ksp_s_atol = 1.0e-10;
  const double ksp_s_dtol = 1.0e5;
  const int ksp_s_maxits = 500; 
  KSPSetTolerances(ksp_s, ksp_s_rtol, ksp_s_atol, ksp_s_dtol, ksp_s_maxits);
  KSPSetType(ksp_s, KSPCG);

  PC pc_s;
  KSPGetPC(ksp_s, &pc_s);
  PCSetType(pc_s, PCJACOBI);

  // Now we add prefix for solver_0, solver_1, and ksp_s
  KSPSetOptionsPrefix(ksp_s, "ls0_inner_");

  KSPSetFromOptions(solver_0->ksp);
  KSPSetFromOptions(solver_1->ksp);
  KSPSetFromOptions(ksp_s);

  // Now we add prefix for pc0, pc1, pc_s
  PCSetOptionsPrefix(pc_s, "pc0_inner_");

  PCSetFromOptions(pc0);
  PCSetFromOptions(pc1);
  PCSetFromOptions(pc_s);

  // print info of solver
  solver_0 -> Info();
  solver_1 -> Info();

  // Copy the RHS
  VecCopy(gAssem_ptr->G_0, G0);
  VecCopy(gAssem_ptr->G_1, G1);

  // Step 1: Solve K_11 yu = G1
  solver_1 -> Solve( G1, sol_1 );

  // Step 2: sol_1 <-- (-1.0)sol_1
  sol_1 -> ScaleValue( -1.0 );

  // Step 3: G0 <-- K_01 sol_1 + G0
  MatMultAdd( gAssem_ptr->K_01, sol_1->solution, G0, G0);

  // Step 4: Solve S yp = G0
  solver_0 -> Solve( G0, sol_0 );

  // Step 5: G1 <- -(K_10 yp - G1)
  VecScale(G1, -1.0);
  MatMultAdd( gAssem_ptr->K_10, sol_0->solution, G1, G1);
  VecScale(G1, -1.0);

  // Step 6: Solve K_11 yu = G1
  solver_1 -> Solve( G1, sol_1 );

  // Clean memory
  MatDestroy(&S);
}


// EOF
