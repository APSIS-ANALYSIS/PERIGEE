#include "PGAssem_LinearPDE_GenAlpha.hpp"

PGAssem_LinearPDE_GenAlpha::PGAssem_LinearPDE_GenAlpha(
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    std::unique_ptr<ALocal_NBC> in_nbc,
    std::unique_ptr<ALocal_EBC> in_ebc,
    std::unique_ptr<IPLocAssem> in_locassem,
    const int &in_nz_estimate )
: locien( std::move(in_locien) ),
  locelem( std::move(in_locelem) ),
  fnode( std::move(in_fnode) ),
  pnode( std::move(in_pnode) ),
  nbc( std::move(in_nbc) ),
  ebc( std::move(in_ebc) ),
  locassem(std::move(in_locassem)),
  num_ebc( ebc->get_num_ebc() ),
  nLocBas( locassem->get_nLocBas() ),
  snLocBas( locassem->get_snLocBas() ),
  dof_mat( locassem->get_dof_mat() ),
  nlgn( pnode->get_nlocghonode() )
{
  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_mat != nbc->get_dof_LID(),
      "PGAssem_NS_FEM::dof_mat != nbc->get_dof_LID(). \n");
  
  const int nlocrow = dof_mat * pnode -> get_nlocalnode();

  // Allocate the sparse matrix K
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat*in_nz_estimate, NULL, dof_mat*in_nz_estimate, NULL, &K);

  // Allocate the vector G
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate();

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);
  
  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_LinearPDE_GenAlpha::~PGAssem_LinearPDE_GenAlpha()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_LinearPDE_GenAlpha::EssBC_KG( const int &field )
{
  const int local_dir = nbc -> get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc->get_LDN(field, ii) * dof_mat + field;

      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(field);
  if(local_sla > 0)
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc->get_LPSN(field, ii) * dof_mat + field;
      const int col = nbc->get_LPMN(field, ii) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_LinearPDE_GenAlpha::EssBC_G( const int &field )
{
  const int local_dir = nbc->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_LinearPDE_GenAlpha::Assem_nonzero_estimate()
{
  const int nElem = locelem->get_nlocalele();
  
  locassem->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = locien -> get_LIEN(ee, ii);
      
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat * ii + mm] = dof_mat * nbc->get_LID( mm, loc_index ) + mm;
    }

    MatSetValues(K, dof_mat * nLocBas, row_index, dof_mat * nLocBas, row_index, locassem->Tangent, ADD_VALUES);
  }
  
  delete [] row_index; row_index = nullptr;
  
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_LinearPDE_GenAlpha::NatBC_G( 
    const double &curr_time, const double &dt )
{
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc -> get_SIEN(ebc_id, ee, LSIEN);

      ebc -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      locassem->Assem_Residual_EBC(ebc_id, curr_time, dt,
          sctrl_x, sctrl_y, sctrl_z);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc -> get_LID(mm, LSIEN[ii]) + mm;
      }

      VecSetValues(G, dof_mat*snLocBas, srow_index, locassem->sur_Residual, ADD_VALUES);
    }
  }

  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

void PGAssem_LinearPDE_GenAlpha::Assem_residual(
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const double &curr_time,
    const double &dt )
{
  const int nElem = locelem->get_nlocalele();

  double * array_a = new double [nlgn * dof_mat];
  double * array_b = new double [nlgn * dof_mat];
  double * local_a = new double [nLocBas * dof_mat];
  double * local_b = new double [nLocBas * dof_mat];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  dot_sol -> GetLocalArray( array_a );
  sol     -> GetLocalArray( array_b );

  for( int ee=0; ee<nElem; ++ee )
  {
    locien -> get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Residual(curr_time, dt, local_a, local_b,
        ectrl_x, ectrl_y, ectrl_z);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc -> get_LID(mm, IEN_e[ii]) + mm;
    }

    VecSetValues(G, dof_mat * nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Resistance type boundary condition
  NatBC_G( curr_time, dt );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_G( ii );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_LinearPDE_GenAlpha::Assem_tangent_residual(
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const double &curr_time,
    const double &dt )
{
  const int nElem = locelem->get_nlocalele();

  double * array_a = new double [nlgn * dof_mat];
  double * array_b = new double [nlgn * dof_mat];
  double * local_a = new double [nLocBas * dof_mat];
  double * local_b = new double [nLocBas * dof_mat];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  dot_sol -> GetLocalArray( array_a );
  sol     -> GetLocalArray( array_b );

  for(int ee=0; ee<nElem; ++ee)
  {
    locien->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
        ectrl_x, ectrl_y, ectrl_z);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat * nbc->get_LID(mm, IEN_e[ii]) + mm;
    }

    MatSetValues(K, dof_mat * nLocBas, row_index, dof_mat * nLocBas, row_index,
        locassem->Tangent, ADD_VALUES);

    VecSetValues(G, dof_mat * nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Natural type boundary condition
  NatBC_G( curr_time, dt );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_LinearPDE_GenAlpha::Assem_mass_residual(
    const PDNSolution * const &sol )
{
  const int nElem = locelem->get_nlocalele();

  double * array_a = new double [nlgn * dof_mat];
  double * local_a = new double [nLocBas * dof_mat];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_a );

  for(int ee=0; ee<nElem; ++ee)
  {
    locien->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Mass_Residual( local_a, ectrl_x, ectrl_y, ectrl_z );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc -> get_LID(mm, IEN_e[ii]) + mm;
    }

    MatSetValues(K, dof_mat * nLocBas, row_index, dof_mat * nLocBas, row_index,
        locassem->Tangent, ADD_VALUES);

    VecSetValues(G, dof_mat * nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Natural type boundary condition
  NatBC_G( 0.0, 0.0 );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

// EOF
