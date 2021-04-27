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

void Matrix_PETSc::gen_ring_inplane_bc( const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dofMat();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_ring_inplane_bc assumes 4 dofs per node.\n");
  
  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      // For non-essential bc nodes, it's diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();
  
  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
    const int dcomp = ring_bc_part -> get_dominant_comp( ii );
    const int row = dnode * dof + dcomp + 1;
    const double bot = -1.0 * ring_bc_part -> get_outvec(ii, dcomp);
    for(int jj=0; jj<3; ++jj)
    {
      if( jj != dcomp )
      {
        const int col = dnode * dof + jj + 1;
        const double top = ring_bc_part -> get_outvec(ii, jj);
        MatSetValue(K, row, col, top/bot, INSERT_VALUES);
      }
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
