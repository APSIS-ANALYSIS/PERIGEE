#include "PLinear_Solver_PETSc.hpp"

PLinear_Solver_PETSc::PLinear_Solver_PETSc()
: rtol( 1.0e-5 ), atol( 1.0e-50 ), dtol( 1.0e50 ), maxits(10000)
{
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  KSPSetFromOptions(ksp);
}

PLinear_Solver_PETSc::PLinear_Solver_PETSc(
    const double &input_rtol, const double &input_atol,
    const double &input_dtol, const int &input_maxits)
: rtol( input_rtol ), atol( input_atol ),
  dtol( input_dtol ), maxits( input_maxits )
{
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  KSPSetFromOptions(ksp);
}

PLinear_Solver_PETSc::PLinear_Solver_PETSc(
    const double &input_rtol, const double &input_atol,
    const double &input_dtol, const int &input_maxits, 
    const char * const &ksp_prefix, const char * const &pc_prefix )
: rtol( input_rtol ), atol( input_atol ),
  dtol( input_dtol ), maxits( input_maxits )
{
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  KSPSetOptionsPrefix( ksp, ksp_prefix );

  PC ksp_pc;
  KSPGetPC( ksp, &ksp_pc );
  PCSetOptionsPrefix( ksp_pc, pc_prefix );
  
  KSPSetFromOptions(ksp);
  PCSetFromOptions(ksp_pc);
}

PLinear_Solver_PETSc::~PLinear_Solver_PETSc()
{
  KSPDestroy(&ksp);
}

void PLinear_Solver_PETSc::Solve( const Vec &G, Vec &out_sol, const bool &isPrint )
{
  KSPSolve(ksp, G, out_sol);

  if( isPrint )
  {
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    PetscReal resnorm;
    KSPGetResidualNorm(ksp, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- KSP: %d, %e", its, resnorm);
  }
}

void PLinear_Solver_PETSc::Solve( const Mat &K, const Vec &G, Vec &out_sol,
   const bool &isPrint )
{
  KSPSetOperators(ksp, K, K);
  KSPSolve(ksp, G, out_sol);

  if( isPrint )
  {
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    PetscReal resnorm;
    KSPGetResidualNorm(ksp, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- KSP: %d, %e", its, resnorm);
  }
}

void PLinear_Solver_PETSc::Solve( const Mat &K, const Vec &G, 
    PDNSolution * const &out_sol, const bool &isPrint )
{
  KSPSetOperators(ksp, K, K);
  KSPSolve(ksp, G, out_sol->solution);

  if( isPrint )
  {
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    PetscReal resnorm;
    KSPGetResidualNorm(ksp, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- KSP: %d, %e", its, resnorm);
  }
  
  // Update solution entries on ghost nodes
  out_sol->GhostUpdate();
}

void PLinear_Solver_PETSc::Solve( const Vec &G, PDNSolution * const &out_sol,
    const bool &isPrint )
{
  KSPSolve(ksp, G, out_sol->solution);

  if( isPrint )
  {
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    PetscReal resnorm;
    KSPGetResidualNorm(ksp, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- KSP: %d, %e ", its, resnorm);
  }
  
  // Update solution entries on ghost nodes
  out_sol->GhostUpdate();
}

void PLinear_Solver_PETSc::Monitor() const
{
#if PETSC_VERSION_LT(3,7,0)
  KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
#elif PETSC_VERSION_LT(3,15,0)
  PetscViewerAndFormat *vf;
  PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
  KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorDefault,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
#else
  PetscViewerAndFormat *vf;
  PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
  KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorResidual,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
#endif
}

//EOF
