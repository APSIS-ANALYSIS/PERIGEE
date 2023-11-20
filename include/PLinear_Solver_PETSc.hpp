#ifndef PLINEAR_SOLVER_PETSC_HPP
#define PLINEAR_SOLVER_PETSC_HPP
// ============================================================================
// PLinear_Solver_PETSc.hpp
//
// Parallel Linear Solver wapper based on PETSc.
//
// Author: Ju Liu
// Date: Dec 8, 2013
// ============================================================================
#include "PDNSolution.hpp"

class PLinear_Solver_PETSc
{
  public:
    KSP ksp;

    // ------------------------------------------------------------------------
    // ! Default KSP object: rtol = 1.0e-5, atol = 1.0e-50, dtol = 1.0e50
    //   maxits = 10000, and the user may reset the parameters by options
    // ------------------------------------------------------------------------
    PLinear_Solver_PETSc();
    
    // ------------------------------------------------------------------------
    // ! Construct KSP with input tolerances and maximum iteration
    // ------------------------------------------------------------------------
    PLinear_Solver_PETSc( const double &input_rtol, const double &input_atol, 
        const double &input_dtol, const int &input_maxits );
    
    // ------------------------------------------------------------------------
    // ! Construct KSP with input tolerances and maximum iteration
    // ------------------------------------------------------------------------
    PLinear_Solver_PETSc( const double &in_rtol, const double &in_atol,
        const double &in_dtol, const int &in_maxits,
        const char * const &ksp_prefix, const char * const &pc_prefix );

    // ------------------------------------------------------------------------
    // Destructor 
    // ------------------------------------------------------------------------
    ~PLinear_Solver_PETSc();

    // ------------------------------------------------------------------------
    // ! Assign a matrix K to the linear solver object
    // ------------------------------------------------------------------------
    void SetOperator(const Mat &K) {KSPSetOperators(ksp, K, K);}

    void SetOperator(const Mat &K, const Mat &P) {KSPSetOperators(ksp, K, P);}

    // ------------------------------------------------------------------------
    // ! Solve a linear problem K out_sol = G
    //   This out_sol is a plain vector, with no ghost entries.
    // ------------------------------------------------------------------------
    void Solve( const Vec &G, Vec &out_sol, const bool &isPrint=true );

    void Solve( const Mat &K, const Vec &G, Vec &out_sol, 
        const bool &isPrint=true );
    
    // ------------------------------------------------------------------------
    // ! Solve a linear problem K out_sol = G
    // ------------------------------------------------------------------------
    void Solve( const Mat &K, const Vec &G, PDNSolution * const &out_sol,
        const bool &isPrint=true );

    // ------------------------------------------------------------------------
    // ! Solve a linear problem with RHS given as G and solution is out_sol
    //   assume the operator has been assigned to KSP
    // ------------------------------------------------------------------------
    void Solve( const Vec &G, PDNSolution * const &out_sol,
        const bool &isPrint=true );

    // ------------------------------------------------------------------------
    // ! Link the solver to a preconditioner context
    // ------------------------------------------------------------------------
    void GetPC( PC *prec ) const {KSPGetPC(ksp, prec);}

    // ------------------------------------------------------------------------
    // ! Get the iteration number
    // ------------------------------------------------------------------------
    int get_ksp_it_num() const
    {
      int it_num; 
      KSPGetIterationNumber(ksp, &it_num); 
      return it_num;
    }

    // ------------------------------------------------------------------------
    // ! Get maximum iteration number for this linear solver
    // ------------------------------------------------------------------------
    int get_ksp_maxits() const
    { 
      int mits; 
#if PETSC_VERSION_LT(3,19,0)
      KSPGetTolerances(ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &mits);
#else
      KSPGetTolerances(ksp, PETSC_NULLPTR, PETSC_NULLPTR, PETSC_NULLPTR, &mits);
#endif
      return mits;
    }

    // ------------------------------------------------------------------------
    // ! Print the ksp info on screen
    // ------------------------------------------------------------------------
    void print_info() const
    {
      SYS_T::print_sep_line();
      KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
      SYS_T::print_sep_line();
    }

    // ------------------------------------------------------------------------
    // ! Monitor the Krylov subspace method behavior
    // ------------------------------------------------------------------------
    void Monitor() const;

  private: 
    // relative, absolute, divergence tolerance
    const PetscReal rtol, atol, dtol;
    
    // maximum number of iterations 
    const PetscInt maxits;
};

#endif
