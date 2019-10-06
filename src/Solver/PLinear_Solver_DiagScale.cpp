#include "PLinear_Solver_DiagScale.hpp"

PLinear_Solver_DiagScale::PLinear_Solver_DiagScale(
    const double &in_rtol, const double &in_atol,
    const double &in_dtol, const int &in_maxits,
    const IPGAssem * const &gAssem_ptr )
: rtol( in_rtol ), atol( in_atol ),
  dtol( in_dtol ), maxits( in_maxits )
{
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  KSPSetFromOptions(ksp);

  // Allocate D
  VecDuplicate( gAssem_ptr->G, &D );
}


PLinear_Solver_DiagScale::PLinear_Solver_DiagScale( 
    const double &in_rtol, const double &in_atol,
    const double &in_dtol, const int &in_maxits,
    const char * const &ksp_prefix, const char * const &pc_prefix,
    const IPGAssem * const &gAssem_ptr )
: rtol( in_rtol ), atol( in_atol ),
  dtol( in_dtol ), maxits( in_maxits )
{
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  KSPSetOptionsPrefix( ksp, ksp_prefix );

  PC ksp_pc;
  KSPGetPC( ksp, &ksp_pc );
  PCSetOptionsPrefix( ksp_pc, pc_prefix );

  KSPSetFromOptions(ksp);
  PCSetFromOptions(ksp_pc);

  // Allocate D
  VecDuplicate( gAssem_ptr->G, &D );
}


PLinear_Solver_DiagScale::~PLinear_Solver_DiagScale()
{
  KSPDestroy(&ksp);
  VecDestroy(&D);
}


void PLinear_Solver_DiagScale::Monitor() const
{
#if PETSC_VERSION_LT(3,7,0)
  KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
#else
  PetscViewerAndFormat *vf;
  PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
  KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorDefault,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
#endif
}


void PLinear_Solver_DiagScale::SymmJacobi_MatVec_Scale( 
    IPGAssem * const &gAssem_ptr )
{
  PETSc_T::GetDiagonal( gAssem_ptr -> K, D ); // obtain diagonal vec

  PETSc_T::MinusSqrtVec( D ); // power to -1/2

  PETSc_T::DiagonalScale( gAssem_ptr-> K, D, D); // symmetric scaling

  PETSc_T::DiagonalScale( D, gAssem_ptr->G ); // Update the G vector
}


void PLinear_Solver_DiagScale::SymmJacobi_Vec_Scale( 
    IPGAssem * const &gAssem_ptr )
{
  PETSc_T::DiagonalScale( D, gAssem_ptr->G );
}


void PLinear_Solver_DiagScale::Solve( const Vec &G, PDNSolution * const &out_sol)
{
#ifdef PETSC_USE_LOG
  PetscLogEvent solver_gmres;
  PetscClassId classid_ls;
  PetscClassIdRegister("petsc_ls", &classid_ls);
  PetscLogEventRegister("petsc solver", classid_ls, &solver_gmres);
  PetscLogEventBegin(solver_gmres,0,0,0,0);
#endif

  KSPSolve(ksp, G, out_sol->solution);

#ifdef PETSC_USE_LOG
  PetscLogEventEnd(solver_gmres,0,0,0,0);
#endif

  KSPGetIterationNumber(ksp, &its);
  KSPGetResidualNorm(ksp, &resnorm);

  PETSc_T::DiagonalScale(D, out_sol->solution);
  out_sol->GhostUpdate();
}


// EOF
