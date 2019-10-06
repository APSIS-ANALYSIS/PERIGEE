#ifndef PLINEAR_SOLVER_DIAGSCALE_SHELL_HPP
#define PLINEAR_SOLVER_DIAGSCALE_SHELL_HPP
// ========================================================
// PLinear_Solver_DiagScale_Shell.hpp
//
// This is a linear solver context derived from 
// PLinear_Solver_DiagScale
// class, which implements a shell pc for the problem.
//
// Author: Ju Liu
// Date: April 7 2018
// ========================================================
#include "PLinear_Solver_DiagScale.hpp"
#include "PEigen_Solver_SLEPc.hpp"

typedef struct {
  // alpha = 1/(rho(B diagF^-1B^tD^-1))
  // gamma = rho(MuF) / 3
  double alpha, gamma;
  
  // D := [diag(Bdiag(F)^(-1)Bt - C)]^-1,
  // Mp := lump(mass_p)^-1
  // Mu := lump(mass_u)^-1
  Vec D, Mp, Mu; 
  
  Mat L, B, C; // L = C Mu B - gamma D
  
  Mat QF; // Mu F

  KSP kspL; // Solver context for L
  
  PEigen_Solver_SLEPc * esolver;
} SchurLSC;

PetscErrorCode SchurPCSetUp( PC, Mat, Mat, Mat, Mat );
PetscErrorCode SchurPCApply( PC, Vec, Vec );

PetscErrorCode SchurPCSetUp_SIMPLE( PC, Mat, Mat, Mat, Mat );
PetscErrorCode SchurPCApply_SIMPLE( PC, Vec, Vec );

PetscErrorCode SchurPCApply_Mass( PC, Vec, Vec );

class PLinear_Solver_DiagScale_Shell : public PLinear_Solver_DiagScale
{
  public:
    PLinear_Solver_DiagScale_Shell( 
        const double &in_rtol, const double &in_atol,
        const double &in_dtol, const int &in_maxits,
        const IPGAssem * const &gAssem_ptr,
        const Mat &Mass_p, const Mat &Mass_u );

    PLinear_Solver_DiagScale_Shell( const double &in_rtol, 
        const double &in_atol, const double &in_dtol, const int &in_maxits,
        const char * const &ksp_prefix, const char * const &pc_prefix,
        const IPGAssem * const &gAssem_ptr,
        const Mat &Mass_p, const Mat &Mass_u );

    virtual ~PLinear_Solver_DiagScale_Shell(); 

    virtual void SetOperator(const Mat &K);

  private:
    SchurLSC * pc_ctx;

};

#endif
