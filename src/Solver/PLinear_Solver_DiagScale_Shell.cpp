#include "PLinear_Solver_DiagScale_Shell.hpp"

PLinear_Solver_DiagScale_Shell::PLinear_Solver_DiagScale_Shell(
    const double &in_rtol, const double &in_atol,
    const double &in_dtol, const int &in_maxits,
    const IPGAssem * const &gAssem_ptr,
    const Mat &Mass_p, const Mat &Mass_u )
: PLinear_Solver_DiagScale( in_rtol, in_atol, in_dtol, 
    in_maxits, gAssem_ptr)
{
  PetscNew( &pc_ctx );

  // Create the PC context
  MatCreateVecs( Mass_p, &(pc_ctx->Mp), NULL);
  MatCreateVecs( Mass_u, &(pc_ctx->Mu), NULL);

  MatGetRowSum( Mass_p, pc_ctx->Mp );
  MatGetRowSum( Mass_u, pc_ctx->Mu );

  VecReciprocal( pc_ctx->Mp );
  VecReciprocal( pc_ctx->Mu );

  VecDuplicate( pc_ctx->Mp, &(pc_ctx->D));
  
  MatDuplicate( Mass_p, MAT_DO_NOT_COPY_VALUES, &(pc_ctx->L) );
  
  MatDuplicate( Mass_u, MAT_DO_NOT_COPY_VALUES, &(pc_ctx->QF) );

  KSPCreate(PETSC_COMM_WORLD, &(pc_ctx -> kspL));
  KSPSetErrorIfNotConverged( pc_ctx->kspL, PETSC_TRUE);
  KSPSetType( pc_ctx->kspL, KSPGMRES );
  // user call slsc_ to access the L solver settings
  KSPSetOptionsPrefix( pc_ctx->kspL, "slsc_");

  // Eigenvalue solver  
  pc_ctx->esolver = new PEigen_Solver_SLEPc();
}


PLinear_Solver_DiagScale_Shell::PLinear_Solver_DiagScale_Shell( 
    const double &in_rtol, const double &in_atol, 
    const double &in_dtol, const int &in_maxits,
    const char * const &ksp_prefix, const char * const &pc_prefix,
    const IPGAssem * const &gAssem_ptr,
    const Mat &Mass_p, const Mat &Mass_u )
: PLinear_Solver_DiagScale( in_rtol, in_atol, in_dtol, 
    in_maxits, ksp_prefix, pc_prefix, gAssem_ptr)
{
  PetscNew( &pc_ctx );

  // Create the PC context
  MatCreateVecs( Mass_p, &(pc_ctx->Mp), NULL);
  MatCreateVecs( Mass_u, &(pc_ctx->Mu), NULL);

  MatGetRowSum( Mass_p, pc_ctx->Mp );
  MatGetRowSum( Mass_u, pc_ctx->Mu );

  VecReciprocal( pc_ctx->Mp );
  VecReciprocal( pc_ctx->Mu );

  VecDuplicate( pc_ctx->Mp, &(pc_ctx->D));
  
  MatDuplicate( Mass_p, MAT_DO_NOT_COPY_VALUES, &(pc_ctx->L) );
  
  MatDuplicate( Mass_u, MAT_DO_NOT_COPY_VALUES, &(pc_ctx->QF) );

  KSPCreate(PETSC_COMM_WORLD, &(pc_ctx -> kspL));
  KSPSetErrorIfNotConverged( pc_ctx->kspL, PETSC_TRUE);
  KSPSetType( pc_ctx->kspL, KSPGMRES );
  // user call slsc_ to access the L solver settings
  KSPSetOptionsPrefix( pc_ctx->kspL, "slsc_");

  // Eigen value solver
  pc_ctx->esolver = new PEigen_Solver_SLEPc();
}


PLinear_Solver_DiagScale_Shell::~PLinear_Solver_DiagScale_Shell()
{
  VecDestroy( &(pc_ctx -> Mp) );
  VecDestroy( &(pc_ctx -> Mu) );
  VecDestroy( &(pc_ctx -> D) );
  MatDestroy( &(pc_ctx -> L) );
  MatDestroy( &(pc_ctx -> QF) );
  KSPDestroy( &(pc_ctx -> kspL) );
  delete pc_ctx->esolver;
  PetscFree( pc_ctx );
}


void PLinear_Solver_DiagScale_Shell::SetOperator(const Mat &K)
{
  KSPSetOperators(ksp, K, K);
  KSPSetUp( ksp );

  PC pc_K;
  KSPGetPC( ksp, &pc_K );

  KSP * sub_ksp;
  PetscInt sub_ksp_num = 2;
  PCFieldSplitGetSubKSP( pc_K, &sub_ksp_num, &sub_ksp );

  PC pc_S;
  KSPGetPC( sub_ksp[1], &pc_S );

  PCSetType( pc_S, PCSHELL );
  PCShellSetContext( pc_S, pc_ctx );
 
  Mat K00, K01, K10, K11; 
  PCFieldSplitGetSchurBlocks( pc_K, &K00, &K01, &K10, &K11 );
  
  SchurPCSetUp( pc_S, K00, K01, K10, K11 );

  PCShellSetApply(pc_S, SchurPCApply);

  PetscFree( sub_ksp );
}


PetscErrorCode SchurPCSetUp( PC spc, Mat Ka, Mat Kb, Mat Kc, Mat Kd )
{
  SchurLSC * ctx;
  PCShellGetContext( spc, (void**)&ctx );

  MatCopy( Ka, ctx->QF, DIFFERENT_NONZERO_PATTERN );  
  MatDiagonalScale( ctx->QF, ctx->Mu, NULL);

  // Assign gamma
  ctx->gamma = ctx->esolver->get_SpectralRad_fast( ctx->QF ) / 3.0;

  // now QF = Mu F Mu
  MatDiagonalScale( ctx->QF, NULL, ctx->Mu );

  // next we use L to store the SIMPLE matrix,
  // L = Kc diag(Ka) Kb - Kd
  MatCreateSchurComplementPmat( Ka, Kb, Kc, Kd,
      MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_REUSE_MATRIX, &(ctx->L) );
  MatScale( ctx->L, -1.0 );

  // Define D
  MatGetDiagonal( ctx->L, ctx->D );

  //PetscReal normD;
  //VecNorm( ctx->D, NORM_INFINITY, &normD );
  
  //Correct gamma
  //ctx->gamma = ctx->gamma / normD;

  // now ctx->D is actually D^-1
  VecReciprocal( ctx->D );

  // L <-- Kc diag(Ka) Kb
  MatAXPY( ctx->L, 1.0, Kd, SUBSET_NONZERO_PATTERN );

  // L <-- Kc diag(Ka) Kb D^-1
  MatDiagonalScale( ctx->L, NULL, ctx->D );

  // define alpha
  ctx->alpha = 1.0 / ctx->esolver->get_SpectralRad_fast( ctx->L );

  // Now ctx->D is alpha D^-1
  VecScale( ctx->D, ctx->alpha );

  // Define L
  MatDestroy( &(ctx->L) );

  Mat MuB;
  MatDuplicate( Kb, MAT_COPY_VALUES, &MuB );
  MatDiagonalScale( MuB, ctx->Mu, NULL );
  MatMatMult( Kc, MuB, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(ctx->L) );
  MatAXPY( ctx->L, -1.0 * ctx->gamma, Kd, DIFFERENT_NONZERO_PATTERN );
  
  KSPSetOperators( ctx->kspL, ctx->L, ctx->L );
  KSPSetFromOptions( ctx->kspL );
 
  ctx->B = Kb;
  ctx->C = Kc;

  MatDestroy(&MuB);
  return 0;
}


PetscErrorCode SchurPCApply( PC spc, Vec sx, Vec sy )
{
  SchurLSC * ctx;

  Vec w, z, u;

  PCShellGetContext( spc, (void**)&ctx );

  MatCreateVecs(ctx->B, NULL, &w);
  MatCreateVecs(ctx->B, NULL, &z);
  MatCreateVecs(ctx->B, &u, NULL);

  KSPSolve( ctx->kspL, sx, sy);
 
  MatMult(ctx->B, sy, w); // w <-- B sy
  MatMult(ctx->QF, w, z); // z <-- Qu^-1 F Qu^-1 w
  MatMult(ctx->C, z, sy); // sy <-- C z

  KSPSolve( ctx->kspL, sy, sy ); // sy <-- L^-1 sy

  VecPointwiseMult(u, sx, ctx->D); // w <-- alpha D^-1 sx

  VecAXPY( sy, -1.0, u );

  VecDestroy(&w);
  VecDestroy(&z);
  VecDestroy(&u);
  return 0;
}


PetscErrorCode SchurPCSetUp_SIMPLE( PC spc, Mat Ka, Mat Kb, Mat Kc, Mat Kd )
{
  SchurLSC * ctx;
  PCShellGetContext( spc, (void**)&ctx );

  // we use L to store the SIMPLE matrix, L = Kd - Kc diag(Ka) Kb
  MatCreateSchurComplementPmat( Ka, Kb, Kc, Kd,
      MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_REUSE_MATRIX, &(ctx->L) );

  KSPSetOperators( ctx->kspL, ctx->L, ctx->L );
  KSPSetFromOptions( ctx->kspL );
 
  return 0;
}


PetscErrorCode SchurPCApply_SIMPLE( PC spc, Vec sx, Vec sy )
{
  SchurLSC * ctx;

  PCShellGetContext( spc, (void**)&ctx );

  KSPSolve( ctx->kspL, sx, sy);
 
  return 0;
}


PetscErrorCode SchurPCApply_Mass( PC spc, Vec sx, Vec sy )
{
  SchurLSC * ctx;

  PCShellGetContext( spc, (void**)&ctx );
  
  VecPointwiseMult(sy, sx, ctx->Mp);
 
  return 0;
}

// EOF
