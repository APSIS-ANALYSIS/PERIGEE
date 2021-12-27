#include "Matrix_PETSc.hpp"

Matrix_PETSc::Matrix_PETSc( const int &loc_row, const int &loc_col,
   const int &dnz, const int &onz )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_col, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  is_set = false;

  MatGetSize(K, &m, &n);

  lm = loc_row;
  ln = loc_col;  
}

Matrix_PETSc::Matrix_PETSc( const int &loc_row, const int &dnz, const int &onz )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_row, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  is_set = false;

  MatGetSize(K, &m, &n);

  lm = loc_row;
  ln = loc_row; 
}

Matrix_PETSc::Matrix_PETSc(const APart_Node * const &pnode_ptr,
    const int &dnz, const int &onz )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  
  lm = pnode_ptr->get_nlocalnode() * pnode_ptr->get_dof();
  ln = lm;
  
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  is_set = false;

  MatGetSize(K, &m, &n);
}

Matrix_PETSc::Matrix_PETSc(const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part, const int &dnz, const int &onz )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  
  lm = pnode_ptr->get_nlocalnode() * bc_part->get_dofMat();
  ln = lm;
  
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  is_set = false;

  MatGetSize(K, &m, &n);
}

Matrix_PETSc::~Matrix_PETSc()
{
  MatDestroy(&K);
}

void Matrix_PETSc::gen_id(const APart_Node * const &pnode_ptr)
{
  if(is_set) Clear();

  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = pnode_ptr->get_dof();
  for(int ii = 0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int index = pnode_ptr->get_node_loc(ii) * dof + jj;
      MatSetValue(K, index, index, 1.0, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc::gen_perm_bc( const APart_Node * const &pnode_ptr,
   const ALocal_NodalBC * const &bc_part )
{
  if(is_set) Clear();
  
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");
  
  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dofMat();
  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc::gen_extractor_for_Dirichlet_nodes( 
    const APart_Node * const &pnode_ptr ,
    const ALocal_NodalBC * const &bc_part )
{
  if(is_set) Clear();

  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = pnode_ptr->get_dof();
  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      if(bc_part->get_LID(jj,ii) == -1)
        MatSetValue(K, row, row, 1.0, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc::MatMultSol( PDNSolution * const &sol ) const
{
  // Check the local dimension of mat and vec
  if(sol->get_nlocal() != ln)
    SYS_T::print_fatal("Error: Matrix_PETSc dimension does not match with the sol vector's dim. \n");

  PDNSolution * temp = new PDNSolution( *sol );

  // temp = K * sol
  MatMult(K, sol->solution, temp->solution);

  // sol = temp
  sol->Copy(*temp);

  // Clean up memory
  delete temp; temp = nullptr;
}

// EOF
