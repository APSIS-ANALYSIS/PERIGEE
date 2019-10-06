#ifndef PEIGEN_SOLVER_SLEPC_HPP
#define PEIGEN_SOLVER_SLEPC_HPP
// ==================================================================
// PEigen_Solver_SLEPc.hpp
//
// Parallel Eigenvalue Solver wrapper based on SLEPc, EPS object.
//
// This class provides an access to the eigenvalue algorithms 
// implemented in the SLEPc library.
//
// Solution method selection --> Sec. 2.4 in SLEPc user manual.
//
// Author: Ju Liu
// Date: April 4 2018
// ==================================================================
#include "Sys_Tools.hpp"
#include "slepceps.h"

class PEigen_Solver_SLEPc
{
  public:
    // The SLEPc object is exposed to the user to ease use.
    EPS eps;

    // Constructor, with tolerance and maximum iteration number
    // Default relative tolerance 1.0e-10, maximum 1000 iterations,
    // and problem is non-Hermitian
    PEigen_Solver_SLEPc( const double &input_rtol = 1.0e-10, 
        const int &input_maxits=1000, const EPSProblemType &ptype=EPS_NHEP );

    // Detructor, Destroy funciton
    virtual ~PEigen_Solver_SLEPc();

    // Print the basic information of the eps object
    virtual void print_info() const;

    // Find the largest eigenvalue in magnitude with the default
    // Krylov-Schur method as the eigen solver method
    double get_SpectralRad( const Mat &K, 
        const EPSType &etype=EPSKRYLOVSCHUR) const;
  
    // Fast calculation of the spectral radius with rough tolerance,
    // less iteration, and do not check the EPSGetConverged 
    double get_SpectralRad_fast( const Mat &K, 
        const double &in_rtol=1.0e-2, const int &in_maxits=100,
        const EPSType &etype=EPSKRYLOVSCHUR ) const;

    // Find the smallest nonzero eigenvalue in magnitude with the default
    // Krylov-Schur method
    double get_SmallestNZ( const Mat &A, const Mat &B,
        const double &in_rtol = 1.0e-6,
        const int &in_maxits = 1000, 
        const EPSType &etype=EPSKRYLOVSCHUR ) const;

    // Print the most recent eigen solver's information
    // Number of iterations,
    // Solution method,
    // Number of eigenvalues and dimension of the subspace, in which
    // nev : number of eigenvalues to compute
    // ncv : the maximum dimension of the subspace to be used by the solver
    // mpd : the maximum dimension allowed for the projected problem
    // Stopping condition
    // Number of converged eigenpairs.
    void print_solve_details() const;
};

#endif
