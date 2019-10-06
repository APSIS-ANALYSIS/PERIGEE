#include "PLinear_Solver_2x2Block_BIPN.hpp"

PLinear_Solver_2x2Block_BIPN::PLinear_Solver_2x2Block_BIPN( const double &rtol0,
    const double &atol0, const double &dtol0, const int &maxit0,
    const double &rtol1, const double &atol1, const double &dtol1,
    const int &maxit1, const double &b_rtol, const int &b_maxits,
    const IPGAssem_2x2Block * const &gAssem_ptr )
: IPLinear_Solver_2x2Block( rtol0, atol0, dtol0, maxit0,
    rtol1, atol1, dtol1, maxit1, gAssem_ptr ),
  bipn_tol( b_rtol ), bipn_maxits( b_maxits )
{
  // Allocate the R_xx vectors
  VecDuplicateVecs( G0, bipn_maxits, &Rcu );
  VecDuplicateVecs( G0, bipn_maxits, &Rcp );
  VecDuplicateVecs( G1, bipn_maxits, &Rmu );
  VecDuplicateVecs( G1, bipn_maxits, &Rmp );

  // Allocate the y0 and y1 vectors
  VecDuplicateVecs( G0, bipn_maxits, &y0 );
  VecDuplicateVecs( G1, bipn_maxits, &y1 );

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

  // Allocate Schur
  MatCreateSchurComplement( gAssem_ptr->K_11, gAssem_ptr->K_11,
      gAssem_ptr->K_10, gAssem_ptr->K_01, gAssem_ptr->K_00, &S);
}


PLinear_Solver_2x2Block_BIPN::~PLinear_Solver_2x2Block_BIPN()
{
  VecDestroyVecs(bipn_maxits, &Rmu);
  VecDestroyVecs(bipn_maxits, &Rmp);
  VecDestroyVecs(bipn_maxits, &Rcu);
  VecDestroyVecs(bipn_maxits, &Rcp);

  VecDestroyVecs(bipn_maxits, &y0);
  VecDestroyVecs(bipn_maxits, &y1);

  delete A; delete Asol;
  delete [] b; b = NULL;
  delete [] alpha; alpha = NULL;
  
  MatDestroy(&S);
}


void PLinear_Solver_2x2Block_BIPN::SetLHS_RHS( 
    IPGAssem_2x2Block const * const &gAssem_ptr )
{
  // Generate Schur complement
  MatCreateSchurComplement( gAssem_ptr->K_11, gAssem_ptr->K_11,
      gAssem_ptr->K_10, gAssem_ptr->K_01, gAssem_ptr->K_00, &S);

  // Assign matrices to solvers
  solver_0 -> SetOperator( S );
  solver_1 -> SetOperator( gAssem_ptr->K_11 );

  // solver_0 is a matrix-free object, cannot do algebraic preconditioner
  PC pc0; solver_0 -> GetPC(&pc0); PCSetType(pc0, PCNONE);

  // Setup the inner solver options in the Schur S
  KSP ksp_s; MatSchurComplementGetKSP(S, &ksp_s);
  KSPSetTolerances(ksp_s, 2.0e-2, 1.0e-20, 1.0e5, 100);

  PC pc_s; KSPGetPC(ksp_s, &pc_s); PCSetType(pc_s, PCJACOBI);

  KSPSetOptionsPrefix(ksp_s, "ls0_inner_");
  PCSetOptionsPrefix(pc_s, "pc0_inner_");

  KSPSetFromOptions(ksp_s);
  PCSetFromOptions(pc_s);

  // set the RHS
  VecCopy(gAssem_ptr->G_0, G0);
  VecCopy(gAssem_ptr->G_1, G1);
}


void PLinear_Solver_2x2Block_BIPN::DestroyOperators()
{
  MatDestroy(&S);
}


void PLinear_Solver_2x2Block_BIPN::Solve_RHS( 
    IPGAssem_2x2Block const * const &gAssem_ptr,
    PDNSolution * const &sol_0, PDNSolution * const &sol_1 )
{
  // Get the initial residual size
  const double G0norm = PETSc_T::Get2Norm(G0);
  const double G1norm = PETSc_T::Get2Norm(G1);
  const double res_e0 = sqrt( G0norm * G0norm + G1norm * G1norm );
  const double res_e0_square = res_e0 * res_e0;
  const double res_e0_tol = res_e0 * bipn_tol;
  double res_ei = res_e0_tol + 1.0;

#ifdef PETSC_USE_LOG
  PetscLogEvent solver_gmres_1, solver_gmres_2, solver_schur, solver_bipn;
  PetscClassId classid_block_ls;
  PetscClassIdRegister("block_ls", &classid_block_ls);
  PetscLogEventRegister("gmres_1", classid_block_ls, &solver_gmres_1);
  PetscLogEventRegister("gmres_2", classid_block_ls, &solver_gmres_2);
  PetscLogEventRegister("solver_schur", classid_block_ls, &solver_schur);
  PetscLogEventRegister("solver_bipn", classid_block_ls, &solver_bipn);
#endif

  // Initialize A and b.
  A -> gen_id();
  for(int ii=0; ii<2*bipn_maxits; ++ii) b[ii] = 0.0;

  // Iteration index
  int ii = 0;

  do
  {
    // 1. GMRES solve Rmp <-- K_11^-1 G1
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solver_gmres_1,0,0,0,0);
    solver_1 -> Solve( G1, Rmp[ii] );
    PetscLogEventEnd(solver_gmres_1,0,0,0,0);
#endif

    // 2. Rmp <- -1.0 * Rmp 
    VecScale( Rmp[ii], -1.0 );

    // 3. Rcp <- K_01 Rmp + G0
    MatMultAdd( gAssem_ptr->K_01, Rmp[ii], G0, Rcp[ii]);

    // 4. Solve S yp = Rcp
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solver_schur,0,0,0,0);
    solver_0 -> Solve( Rcp[ii], y0[ii] );
    PetscLogEventEnd(solver_schur,0,0,0,0);
#endif

    // 5. Rmp := K_10 yp
    MatMult( gAssem_ptr->K_10, y0[ii], Rmp[ii] ); 

    // 6. Rmu = Rm - Rmp
    VecWAXPY(Rmu[ii], -1.0, Rmp[ii], G1);

    // 7. Solve K11 yu = Rmu
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solver_gmres_2,0,0,0,0);
    solver_1 -> Solve( Rmu[ii], y1[ii] );
    PetscLogEventEnd(solver_gmres_2,0,0,0,0);
#endif

    // 8. Calculate Rmu, Rcp, Rcu
    MatMult( gAssem_ptr->K_11, y1[ii], Rmu[ii] );
    MatMult( gAssem_ptr->K_00, y0[ii], Rcp[ii] );
    MatMult( gAssem_ptr->K_01, y1[ii], Rcu[ii] );

    // 9. Assemble the A matrix
    // aa and bb are two temporary double for holding VecDot result
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solver_bipn,0,0,0,0);
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
    PetscLogEventEnd(solver_bipn,0,0,0,0);
#endif

    // 13. Calculate the error at this iteration
    res_ei = 0.0;
    for(int jj=0; jj<=ii; ++jj)
      res_ei += alpha[jj]*b[jj] + alpha[jj+bipn_maxits]*b[jj+bipn_maxits];

    // In the original paper, the algorithm is res_e0_sqaure - res_ei.
    // When the difference is tiny, it can be very small negative values.
    // So I use std::abs to make this step more accurate.
    res_ei = sqrt( std::abs(res_e0_square - res_ei) );

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
}


// EOF
