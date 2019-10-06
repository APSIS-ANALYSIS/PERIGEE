#ifndef PLINEAR_SOLVER_DIAGSCALE_HPP
#define PLINEAR_SOLVER_DIAGSCALE_HPP
// ==================================================================
// PLinear_Solver_DiagScale.hpp
//
// This is a linear solver class that implements explicitly a 
// diagonal scaled solver.
//
// A Vec object W is used as the diag(K)^-0.5.
//
// For the original matrix problem Kx = b,
// a new matrix problem DKDy = Db is solved first and then
//                         x = Dy.
//
// This class is used as a fair comparison for the GMRES with BIPN.
//
// Author: Ju Liu
// Date: Mar. 8 2018
// ==================================================================
#include "PETSc_Tools.hpp"
#include "IPGAssem.hpp"

class PLinear_Solver_DiagScale
{
  public:
    KSP ksp;

    PLinear_Solver_DiagScale( const double &in_rtol, const double &in_atol,
        const double &in_dtol, const int &in_maxits,
        const IPGAssem * const &gAssem_ptr );

    PLinear_Solver_DiagScale( const double &in_rtol, const double &in_atol,
        const double &in_dtol, const int &in_maxits,
        const char * const &ksp_prefix, const char * const &pc_prefix,
        const IPGAssem * const &gAssem_ptr );

    virtual ~PLinear_Solver_DiagScale();

    // --------------------------------------------------------------
    // D <-- Diag(K)^-0.5,
    // K <-- D K D,
    // G <-- D G.
    // Note: This will modify the gAssem_ptr's K and G objects.
    // --------------------------------------------------------------
    virtual void SymmJacobi_MatVec_Scale( IPGAssem * const &gAssem_ptr );

    // --------------------------------------------------------------
    // Scale the RHS vector by G <-- D G.
    // One should make sure that the D has been obtained by the 
    // SymmJacobi_MatVec_Scale function.
    // --------------------------------------------------------------
    virtual void SymmJacobi_Vec_Scale( IPGAssem * const &gAssem_ptr );


    virtual void SetOperator(const Mat &K) {KSPSetOperators(ksp, K, K);}

    // --------------------------------------------------------------
    // Link the solver to a preconditioner context
    // --------------------------------------------------------------
    virtual void GetPC( PC *prec ) const {KSPGetPC(ksp, prec);}

    // --------------------------------------------------------------
    // Solve for K y = G, and then update solution as x = Dy.
    // --------------------------------------------------------------
    virtual void Solve( const Vec &G, PDNSolution * const &out_sol);


    virtual int get_ksp_it_num() const
    {int it_num; KSPGetIterationNumber(ksp, &it_num); return it_num;}


    virtual void Print_ConvInfo() const
    {
      PetscReal rnorm;
      KSPGetResidualNorm(ksp, &rnorm);
      PetscPrintf(PETSC_COMM_WORLD, " KSP Residual: %e ", rnorm);
    }

    virtual void Monitor() const;

    virtual void Info() const
    {
      SYS_T::commPrint("----------------------------------------------------------- \n");
      KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
      SYS_T::commPrint("----------------------------------------------------------- \n");
    }

  private:
    const PetscReal rtol, atol, dtol;

    const PetscInt maxits;

    PetscInt its;

    PetscReal resnorm;

    // This is the additional data that stores the diag(K)^-0.5
    Vec D;
};

#endif
