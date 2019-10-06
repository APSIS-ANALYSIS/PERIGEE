#include "PLinear_Solver_BIPN_2x2Block.hpp"

PLinear_Solver_BIPN_2x2Block::PLinear_Solver_BIPN_2x2Block(
    const double &rtol0, const double &atol0,
    const double &dtol0, const int &maxit0,
    const double &rtol1, const double &atol1,
    const double &dtol1, const int &maxit1,
    const double &b_rtol, const int &b_maxits,
    const IPGAssem_2x2Block * const &gAssem_ptr )
: bipn_tol( b_rtol ), bipn_maxits( b_maxits )
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

  // Allocate the R_xx vectors
  VecDuplicateVecs( G0, bipn_maxits, &Rcu );
  VecDuplicateVecs( G0, bipn_maxits, &Rcp );
  VecDuplicateVecs( G1, bipn_maxits, &Rmu );
  VecDuplicateVecs( G1, bipn_maxits, &Rmp );

  // Allocate the y0 and y1 vectors
  VecDuplicateVecs( G0, bipn_maxits, &y0 );
  VecDuplicateVecs( G1, bipn_maxits, &y1 );

  // Generate an identity matrix with same layout as K_11
  //PETSc_T::MatCreateId( Id, PETSc_T::GetNumLocalRow( gAssem_ptr->K_11 ) );
  
  // Allocate the normal equation and RHS
  A = new MATH_T::Matrix_SymPos_dense( 2*bipn_maxits );
  A -> gen_id();
  Asol = new MATH_T::Matrix_SymPos_dense( 2*bipn_maxits );
  Asol -> gen_id();
  
  b = new double [2*bipn_maxits];
  alpha = new double [2*bipn_maxits];
  
  for(int ii=0; ii<2*bipn_maxits; ++ii) 
  {
    b[ii] = 0.0;
    alpha[ii] = 0.0;
  }
}


PLinear_Solver_BIPN_2x2Block::~PLinear_Solver_BIPN_2x2Block()
{
  delete solver_0; delete solver_1;
  VecDestroy(&G0); VecDestroy(&G1);
  VecDestroy(&D0); VecDestroy(&D1);
  
  VecDestroyVecs(bipn_maxits, &Rmu); 
  VecDestroyVecs(bipn_maxits, &Rmp);
  VecDestroyVecs(bipn_maxits, &Rcu); 
  VecDestroyVecs(bipn_maxits, &Rcp);
  
  VecDestroyVecs(bipn_maxits, &y0);
  VecDestroyVecs(bipn_maxits, &y1);
  
  //MatDestroy(&Id);
  delete A; delete Asol;
  delete [] b; b = NULL;
  delete [] alpha; alpha = NULL;
}


void PLinear_Solver_BIPN_2x2Block::SymmJacobi_Scale(
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


void PLinear_Solver_BIPN_2x2Block::SymmJacobi_Restore(
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


void PLinear_Solver_BIPN_2x2Block::Solve_BIPN(
    IPGAssem_2x2Block * const &gAssem_ptr,
    PDNSolution * const &sol_0, PDNSolution * const &sol_1 )
{
  // Apply symmetric Jacobi preconditioner
  SymmJacobi_Scale(gAssem_ptr); 

  // Copy the RHS
  VecCopy(gAssem_ptr->G_0, G0);
  VecCopy(gAssem_ptr->G_1, G1);

  // Create the SIMPLE Schur complement
  // S := K_00 - K_01 * inv(Id) * K_10
  MatCreateSchurComplement( gAssem_ptr->K_11, gAssem_ptr->K_11,
      gAssem_ptr->K_10, gAssem_ptr->K_01, gAssem_ptr->K_00, &S);

  // Assign operators
  solver_0 -> SetOperator( S );
  solver_1 -> SetOperator( gAssem_ptr->K_11 );

  // solver_0 for this one is a matrix-free fashion, cannot apply
  // algebraic preconditioners.
  PC pc0; solver_0 -> GetPC(&pc0); PCSetType(pc0, PCNONE);

  // Setup the inversion for the identity matrix 
  KSP ksp_s; MatSchurComplementGetKSP(S, &ksp_s);
  KSPSetTolerances(ksp_s, 2.0e-2, 1.0e-20, 1.0e5, 100); // allow at most 2 its

  PC pc_s; KSPGetPC(ksp_s, &pc_s); PCSetType(pc_s, PCJACOBI);

  KSPSetOptionsPrefix(ksp_s, "ls0_inner_");
  PCSetOptionsPrefix(pc_s, "pc0_inner_");

  // Allow the users to change inner solvers
  KSPSetFromOptions(ksp_s);
  PCSetFromOptions(pc_s);

  // Get the initial residual size
  const double G0norm = PETSc_T::Get2Norm(G0);
  const double G1norm = PETSc_T::Get2Norm(G1);
  const double res_e0 = sqrt( G0norm * G0norm + G1norm * G1norm );
  const double res_e0_square = res_e0 * res_e0;
  const double res_e0_tol = res_e0 * bipn_tol;
  double res_ei = res_e0_tol + 1.0;

  // Initialize A and b.
  A -> gen_id();
  for(int ii=0; ii<2*bipn_maxits; ++ii) b[ii] = 0.0;

  // Iteration index
  int ii = 0;

  do
  {
    // 1. GMRES solve Rmp <-- K_11^-1 G1
    solver_1 -> Solve( G1, Rmp[ii] );
    
    // 2. Rmp <- -1.0 * Rmp 
    VecScale( Rmp[ii], -1.0 );

    // 3. Rcp <- K_01 Rmp + G0
    MatMultAdd( gAssem_ptr->K_01, Rmp[ii], G0, Rcp[ii]);

    // 4. Solve S yp = Rcp
    solver_0 -> Solve( Rcp[ii], y0[ii] );

    // 5. Rmp := K_10 yp
    MatMult( gAssem_ptr->K_10, y0[ii], Rmp[ii] ); 

    // 6. Rmu = Rm - Rmp
    VecWAXPY(Rmu[ii], -1.0, Rmp[ii], G1);

    // 7. Solve K11 yu = Rmu
    solver_1 -> Solve( Rmu[ii], y1[ii] );

    // 8. Calculate Rmu, Rcp, Rcu
    MatMult( gAssem_ptr->K_11, y1[ii], Rmu[ii] );
    MatMult( gAssem_ptr->K_00, y0[ii], Rcp[ii] );
    MatMult( gAssem_ptr->K_01, y1[ii], Rcu[ii] );

    // 9. Assemble the A matrix
    // aa and bb are two temporary double for holding VecDot result
    double aa, bb;
    for(int jj=0; jj<=ii; ++jj)
    {
      VecDot(Rmu[ii], Rmu[jj], &aa); VecDot(Rcu[ii], Rcu[jj], &bb);
      A->set_value(ii, jj, aa + bb);
      A->set_value(jj, ii, aa + bb);

      VecDot(Rmu[ii], Rmp[jj], &aa); VecDot(Rcu[ii], Rcp[jj], &bb);
      A->set_value(ii, jj+bipn_maxits, aa+bb);
      A->set_value(jj+bipn_maxits, ii, aa+bb);

      VecDot(Rmp[ii], Rmu[jj], &aa); VecDot(Rcp[ii], Rcu[jj], &bb);
      A->set_value(ii+bipn_maxits, jj, aa+bb);
      A->set_value(jj, ii+bipn_maxits, aa+bb);

      VecDot(Rmp[ii], Rmp[jj], &aa); VecDot(Rcp[ii], Rcp[jj], &bb);
      A->set_value(ii+bipn_maxits, jj+bipn_maxits, aa+bb);
      A->set_value(jj+bipn_maxits, ii+bipn_maxits, aa+bb);
    }
    
    // 10. Update the residual
    VecCopy(gAssem_ptr->G_0, G0);
    VecCopy(gAssem_ptr->G_1, G1);
   
    // 11. Get b vector using the true RHS vectors from global assem ptr 
    VecDot( Rmu[ii], G1, &aa ); VecDot( Rcu[ii], G0, &bb );
    b[ii] = aa + bb;

    VecDot( Rmp[ii], G1, &aa ); VecDot( Rcp[ii], G0, &bb );
    b[bipn_maxits + ii] = aa + bb;

    // 12. Solve for A alpha = b
    Asol -> copy( A );
    Asol -> LDLt_fac();
    Asol -> LDLt_solve( b, alpha );

    // 13. Calculate the error at this iteration
    res_ei = 0.0;
    for(int jj=0; jj<=ii; ++jj)
      res_ei += alpha[jj]*b[jj] + alpha[jj+bipn_maxits]*b[jj+bipn_maxits];
    
    res_ei = sqrt( res_e0_square - res_ei );

    // 14. Recall from step 10, G0 G1 are the true RHS
    //     Now use G0 and G1 to get Rm and Rc for the first stage
    //     Schur complement solve. 
    for(int jj=0; jj<=ii; ++jj)
    {
      VecAXPY( G0, -1.0 * alpha[jj], Rcu[jj] );
      VecAXPY( G0, -1.0 * alpha[jj+bipn_maxits], Rcp[jj] );
      VecAXPY( G1, -1.0 * alpha[jj], Rmu[jj] );
      VecAXPY( G1, -1.0 * alpha[jj+bipn_maxits], Rmp[jj] );
    } 

    // 15. Update the interation number
    ii += 1;
    PetscPrintf(PETSC_COMM_WORLD, " %e \n", res_ei);
  }while( ii<bipn_maxits && res_e0_tol < res_ei );

  // Correct the iteration index if reach the maximum value 
  if(ii == bipn_maxits) ii = bipn_maxits-1; 

  PetscPrintf(PETSC_COMM_WORLD, "  --- BIPN: %d, %e, %e \n", 
      ii, res_ei/res_e0, res_ei);

  // Now Give the acutal solution
  VecSet( sol_0 -> solution , 0.0 );
  VecSet( sol_1 -> solution , 0.0 );

  for(int jj=0; jj<=ii; ++jj)
  {
    VecAXPY( sol_0->solution, alpha[jj+bipn_maxits], y0[jj]);
    VecAXPY( sol_1->solution, alpha[jj], y1[jj]);
  }

  // Apply the Symmetric Jacobi preconditioner
  PETSc_T::DiagonalScale(D0, sol_0->solution);
  PETSc_T::DiagonalScale(D1, sol_1->solution);

  // Update ghost values
  sol_0 -> GhostUpdate();
  sol_1 -> GhostUpdate();

  // Destroy the SIMPLE Schur complement
  MatDestroy(&S);

  // Restore the original matrices and vectors in gAssem_ptr
  SymmJacobi_Restore(gAssem_ptr);
}

// EOF
