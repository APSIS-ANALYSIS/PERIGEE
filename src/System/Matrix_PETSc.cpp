#include "Matrix_PETSc.hpp"

Matrix_PETSc::Matrix_PETSc( const int &loc_row, const int &loc_col,
   const int &dnz, const int &onz ) 
: lm( loc_row ), ln( loc_col ), is_set(false)
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");

#if PETSC_VERSION_LT(3,19,0)
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_col, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);
#else
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_col, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULLPTR, onz, PETSC_NULLPTR, &K);
#endif

  MatGetSize(K, &gm, &gn);
}

Matrix_PETSc::Matrix_PETSc( const int &loc_row, const int &dnz, const int &onz )
: lm( loc_row ), ln( loc_row ), is_set(false)
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");

#if PETSC_VERSION_LT(3,19,0)
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_row, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);
#else
  MatCreateAIJ(PETSC_COMM_WORLD, loc_row, loc_row, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULLPTR, onz, PETSC_NULLPTR, &K);
#endif

  MatGetSize(K, &gm, &gn);
}

Matrix_PETSc::Matrix_PETSc(const APart_Node * const &pnode_ptr,
    const int &dnz, const int &onz ) : is_set( false )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  
  lm = pnode_ptr->get_nlocalnode() * pnode_ptr->get_dof();
  ln = lm;
  
#if PETSC_VERSION_LT(3,19,0)
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);
#else
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULLPTR, onz, PETSC_NULLPTR, &K);
#endif

  MatGetSize(K, &gm, &gn);
}

Matrix_PETSc::Matrix_PETSc(const APart_Node * const &pnode_ptr,
    const ALocal_NBC * const &bc_part, const int &dnz, const int &onz )
: is_set( false )
{
  SYS_T::commPrint("===> PETSc: MatCreateAIJ called. \n");
  
  lm = pnode_ptr->get_nlocalnode() * bc_part->get_dof_LID();
  ln = lm;
  
#if PETSC_VERSION_LT(3,19,0)
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULL, onz, PETSC_NULL, &K);
#else
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DECIDE, PETSC_DECIDE, 
      dnz, PETSC_NULLPTR, onz, PETSC_NULLPTR, &K);
#endif

  MatGetSize(K, &gm, &gn);
}

Matrix_PETSc::~Matrix_PETSc()
{
  MatDestroy(&K);
}

void Matrix_PETSc::gen_id(const APart_Node * const &pnode_ptr)
{
  if(is_set) Clear();

  SYS_T::print_fatal_if(gm != gn, "Error: This is not a square matrix. \n");

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

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

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
   const ALocal_NBC * const &bc_part )
{
  if(is_set) Clear();
  
  SYS_T::print_fatal_if(gm != gn, "Error: This is not a square matrix. \n");
  
  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();
  
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

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

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

void Matrix_PETSc::gen_perm_bc( const std::vector<APart_Node *> &pnode_list,
    const std::vector<ALocal_NBC *> &bc_part_list, 
    const std::vector<int> &start_idx )
{
  if(is_set) Clear();

  SYS_T::print_fatal_if( pnode_list.size() != bc_part_list.size(), "Error: the input apart_node and alocal_nodalbc should have the same length. \n");
  SYS_T::print_fatal_if(gm != gn, "Error: This is not a square matrix. \n");

  const int nfield = VEC_T::get_size(pnode_list);
  
  for(int ff=0; ff<nfield; ++ff)
  {
    const int dof   = pnode_list[ff] -> get_dof();
    const int nnode = pnode_list[ff] -> get_nlocalnode();
    
    for(int ii=0; ii<nnode; ++ii)
    {
      for(int jj=0; jj<dof; ++jj)
      {
        const int row = start_idx[ff] + ii * dof + jj;
        const int col = bc_part_list[ff] -> get_LID(jj, ii);
        MatSetValue(K, row, col, 1.0, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  
  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

  for(int ff=0; ff<nfield; ++ff)
  {
    const int dof   = pnode_list[ff] -> get_dof();
    const int nnode = pnode_list[ff] -> get_nlocalnode();
    for(int ii=0; ii<nnode; ++ii)
    {
      for(int jj=0; jj<dof; ++jj)
      {
        const int row = start_idx[ff] + ii * dof + jj;
        const int col = bc_part_list[ff] -> get_LID(jj, ii);
        MatSetValue(K, row, col, 1.0, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc::gen_extractor_for_Dirichlet_nodes( 
    const APart_Node * const &pnode_ptr ,
    const ALocal_NBC * const &bc_part )
{
  if(is_set) Clear();

  SYS_T::print_fatal_if(gm != gn, "Error: This is not a square matrix. \n");

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

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, lm, ln, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

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
  SYS_T::print_fatal_if( sol->get_nlocal() != ln, "Error: Matrix_PETSc dimension does not match with the sol vector's dim. \n" );

  PDNSolution * temp = new PDNSolution( *sol );

  // temp = K * sol
  MatMult(K, sol->solution, temp->solution);

  // sol = temp
  sol->Copy(*temp);

  // Clean up memory
  delete temp; temp = nullptr;
}

void Matrix_PETSc::MatMultSol( Vec &sol ) const
{
  Vec temp;
  VecDuplicate(sol, &temp);

  MatMult(K, sol, temp);

  VecCopy(temp, sol);

  VecDestroy(&temp);
}

// EOF
