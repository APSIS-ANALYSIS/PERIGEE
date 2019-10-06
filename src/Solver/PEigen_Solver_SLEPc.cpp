#include "PEigen_Solver_SLEPc.hpp"

PEigen_Solver_SLEPc::PEigen_Solver_SLEPc( const double &input_rtol, 
    const int &input_maxits, const EPSProblemType &ptype )
{
  EPSCreate(PETSC_COMM_WORLD, &eps);

  // Set default tolerance and max number of iterations
  EPSSetTolerances( eps, input_rtol, input_maxits );

  // Default the problem to be non-Hermitian
  EPSSetProblemType(eps, ptype);
}


PEigen_Solver_SLEPc::~PEigen_Solver_SLEPc()
{
  EPSDestroy(&eps);
}


void PEigen_Solver_SLEPc::print_info() const
{
  EPSView( eps, PETSC_VIEWER_STDOUT_WORLD );
}


double PEigen_Solver_SLEPc::get_SpectralRad( const Mat &K, 
    const EPSType &etype ) const
{
  // Assign the matrix
  EPSSetOperators(eps, K, NULL);

  // Setup the eigen solver
  EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
  EPSSetType(eps, etype);
  EPSSetDimensions( eps, 1, 20, 20 );

  // Sove the eigenvalue problem
  EPSSolve(eps);
  
  PetscInt nconv;
  EPSGetConverged(eps,&nconv);

  SYS_T::print_fatal_if(nconv == 0, "get_spectrum_radius does not converge.\n");

  // Get the real and imaginary part of the eigenvalue
  PetscScalar kr, ki;
  EPSGetEigenvalue(eps, 0, &kr, &ki);

  return std::sqrt(kr*kr + ki*ki);
}


double PEigen_Solver_SLEPc::get_SpectralRad_fast( const Mat &K, 
    const double &in_rtol, const int &in_maxits, const EPSType &etype ) const
{
  // Assign the matrix
  EPSSetOperators(eps, K, NULL);

  // Setup the eigen solver
  EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
  EPSSetType(eps, etype);
  EPSSetDimensions( eps, 1, 20, 20 );
  EPSSetTolerances( eps, in_rtol, in_maxits );

  // Sove the eigenvalue problem
  EPSSolve(eps);
  
  // Get the real and imaginary part of the eigenvalue
  PetscScalar kr, ki;
  EPSGetEigenvalue(eps, 0, &kr, &ki);

  return std::sqrt(kr*kr + ki*ki);
}


double PEigen_Solver_SLEPc::get_SmallestNZ( const Mat &A, const Mat &B, 
    const double &in_rtol, const int &in_maxits, const EPSType &etype ) const
{
  // Assign the matrix
  EPSSetOperators(eps, A, B);
  
  // Setup the eigen solver
  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
  EPSSetType(eps, etype);
  const int user_nev = 10;
  const int user_ncv = 20;
  const int user_mpd = 20;
  EPSSetDimensions( eps, user_nev, user_ncv, user_mpd );
  EPSSetTolerances( eps, in_rtol, in_maxits );
  
  // Sove the eigenvalue problem
  EPSSolve(eps);
  
  // Get the real and imaginary part of the eigenvalue
  PetscScalar kr, ki, val;
  for(int ii=0; ii<user_nev; ++ii)
  {
    EPSGetEigenvalue(eps, ii, &kr, &ki);
    val = std::sqrt(kr*kr + ki*ki);
    if( std::abs(val) > 0.0 ) break;
  }

  return val;
}


void PEigen_Solver_SLEPc::print_solve_details() const
{
  PetscPrintf(PETSC_COMM_WORLD," \n ------ SLEPc solver info ------ \n");
  PetscInt its;
  EPSGetIterationNumber(eps,&its);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);

  EPSType type; 
  EPSGetType(eps,&type);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);

  PetscInt nev, ncv, mpd; 
  EPSGetDimensions(eps, &nev, &ncv, &mpd);
  PetscPrintf(PETSC_COMM_WORLD,
      " Number of requested eigenvalues: nev = %d, ncv = %d, mpd = %d\n",
      nev, ncv, mpd);

  PetscReal tol; PetscInt maxit; 
  EPSGetTolerances(eps, &tol, &maxit);
  PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",
      (double)tol, maxit);

  PetscInt nconv;
  EPSGetConverged(eps,&nconv);
  PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %d\n",nconv);
  PetscPrintf(PETSC_COMM_WORLD," ------------------------------- \n");
}

// EOF
