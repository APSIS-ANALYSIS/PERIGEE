#include "IPLinear_Solver_2x2Block.hpp"

IPLinear_Solver_2x2Block::IPLinear_Solver_2x2Block(
    const double &rtol0, const double &atol0,
    const double &dtol0, const int &maxit0,
    const double &rtol1, const double &atol1,
    const double &dtol1, const int &maxit1,
    const IPGAssem_2x2Block * const &gAssem_ptr )
{
  // Setup solvers
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


IPLinear_Solver_2x2Block::~IPLinear_Solver_2x2Block()
{
  delete solver_0; solver_0 = NULL;
  delete solver_1; solver_1 = NULL;
  
  VecDestroy(&G0); VecDestroy(&G1);
  VecDestroy(&D0); VecDestroy(&D1);
}


void IPLinear_Solver_2x2Block::SymmJacobi_MatVec_Scale( 
    IPGAssem_2x2Block * const &gAssem_ptr )
{
  // Obtain Diagonal entries
  PETSc_T::GetDiagonal(gAssem_ptr->K_00, D0);
  PETSc_T::GetDiagonal(gAssem_ptr->K_11, D1);

  // Power to -1/2
  PETSc_T::MinusSqrtVec( D0 );
  PETSc_T::MinusSqrtVec( D1 );

  // Update the K matrices in gAssem_ptr
  PETSc_T::DiagonalScale(gAssem_ptr->K_00, D0, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_01, D0, D1);
  PETSc_T::DiagonalScale(gAssem_ptr->K_10, D1, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_11, D1, D1);

  // Update the G vectors in gAssem_ptr
  PETSc_T::DiagonalScale( D0, gAssem_ptr->G_0 );
  PETSc_T::DiagonalScale( D1, gAssem_ptr->G_1 );
}


void IPLinear_Solver_2x2Block::SymmJacobi_MatVec_Restore( 
    IPGAssem_2x2Block * const &gAssem_ptr )
{
  // Power to -1
  PETSc_T::InvAbsVec( D0 );
  PETSc_T::InvAbsVec( D1 ); 

  // Update K matrices
  PETSc_T::DiagonalScale(gAssem_ptr->K_00, D0, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_01, D0, D1);
  PETSc_T::DiagonalScale(gAssem_ptr->K_10, D1, D0);
  PETSc_T::DiagonalScale(gAssem_ptr->K_11, D1, D1);

  // Update G vectors
  PETSc_T::DiagonalScale( D0, gAssem_ptr->G_0 );
  PETSc_T::DiagonalScale( D1, gAssem_ptr->G_1 );
}


// EOF
