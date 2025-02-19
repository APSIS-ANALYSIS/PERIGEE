#include "PGAssem_NS_FEM.hpp"

PGAssem_NS_FEM::PGAssem_NS_FEM(
    // IPLocAssem * const &locassem_ptr,
    // FEAElement * const &elements,
    // const IQuadPts * const &quads,
    // const AGlobal_Mesh_Info * const &agmi_ptr,
    // const ALocal_Elem * const &alelem_ptr,
    // const ALocal_IEN * const &aien_ptr,
    // const APart_Node * const &pnode_ptr,
    // const ALocal_NBC * const &part_nbc,
    const ALocal_EBC * const &part_ebc,
    const IGenBC * const &gbc,
    // const int &in_nz_estimate )
// : nLocBas( agmi_ptr->get_nLocBas() ),
//   dof_sol( pnode_ptr->get_dof() ),
//   dof_mat( locassem_ptr->get_dof_mat() ),
//   num_ebc( part_ebc->get_num_ebc() ),
//   nlgn( pnode_ptr->get_nlocghonode() ),
//   snLocBas( 0 ) 
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    // std::unique_ptr<ALocal_InflowBC> in_infbc,
    std::unique_ptr<ALocal_NBC> in_nbc,
    // std::unique_ptr<ALocal_EBC> in_ebc,
    // std::unique_ptr<IGenBC> in_gbc,
    std::unique_ptr<ALocal_WeakBC> in_wbc,
    std::unique_ptr<IPLocAssem> in_locassem,    
    const int &in_nz_estimate )
: locien( std::move(in_locien) ),
  locelem( std::move(in_locelem) ),
  fnode( std::move(in_fnode) ),
  pnode( std::move(in_pnode) ),
  // infbc( std::move(in_infbc) ),
  nbc( std::move(in_nbc) ),
  // ebc( std::move(in_ebc) ),
  // gbc( std::move(in_gbc) ),
  wbc( std::move(in_wbc) ),
  locassem(std::move(in_locassem)),
  nLocBas( locassem->get_nLocBas() ),
  snLocBas( locassem->get_snLocBas() ),
  dof_sol( pnode->get_dof() ),
  dof_mat( locassem->get_dof_mat() ),
  num_ebc( ebc->get_num_ebc() ),
  nlgn( pnode->get_nlocghonode() )
{
  // Make sure the data structure is compatible
  // SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),

  // SYS_T::print_fatal_if(dof_mat != part_nbc->get_dof_LID(),
  //     "PGAssem_NS_FEM::dof_mat != part_nbc->get_dof_LID(). \n");
  //     "PGAssem_NS_FEM::dof_sol != locassem_ptr->get_dof(). \n");
  SYS_T::print_fatal_if(dof_sol != locassem->get_dof(),
      "PGAssem_NS_FEM::dof_sol != locassem->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != nbc->get_dof_LID(),
      "PGAssem_NS_FEM::dof_mat != nbc->get_dof_LID(). \n");

  // Make sure that the surface element's number of local basis are 
  // the same. This is an assumption in this assembly routine.
  // if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_NS_FEM, snLocBas has to be uniform. \n");
  }

  const int nlocrow = dof_mat * pnode->get_nlocalnode();

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

  // Assem_nonzero_estimate( alelem_ptr, locassem_ptr, 
  //     elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ebc, gbc );
  // Assem_nonzero_estimate();
  Assem_nonzero_estimate( part_ebc, gbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation
 
  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_NS_FEM::~PGAssem_NS_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

// void PGAssem_NS_FEM::EssBC_KG(
//     const ALocal_NBC * const &nbc_part, const int &field )
void PGAssem_NS_FEM::EssBC_KG( const int &field )
{
  // const int local_dir = nbc_part->get_Num_LD(field);
  const int local_dir = nbc->get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      // const int row = nbc_part->get_LDN(field, i) * dof_mat + field;
      const int row = nbc->get_LDN(field, i) * dof_mat + field;
      
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  // const int local_sla = nbc_part->get_Num_LPS(field);
  const int local_sla = nbc->get_Num_LPS(field);
  if(local_sla > 0)
  {
    for(int i=0; i<local_sla; ++i)
    {
    //   const int row = nbc_part->get_LPSN(field, i) * dof_mat + field;
    //   const int col = nbc_part->get_LPMN(field, i) * dof_mat + field;
      const int row = nbc->get_LPSN(field, i) * dof_mat + field;
      const int col = nbc->get_LPMN(field, i) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

// void PGAssem_NS_FEM::EssBC_G( const ALocal_NBC * const &nbc_part, 
//     const int &field )
void PGAssem_NS_FEM::EssBC_G( const int &field )
{
  // const int local_dir = nbc_part->get_Num_LD(field);
  const int local_dir = nbc->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      // const int row = nbc_part->get_LDN(field, ii) * dof_mat + field;
      const int row = nbc->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  // const int local_sla = nbc_part->get_Num_LPS(field);
  const int local_sla = nbc->get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      // const int row = nbc_part->get_LPSN(field, ii) * dof_mat + field;
      const int row = nbc->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

// void PGAssem_NS_FEM::Assem_nonzero_estimate(
//     const ALocal_Elem * const &alelem_ptr,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &elements,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const APart_Node * const &node_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const IGenBC * const &gbc )
void PGAssem_NS_FEM::Assem_nonzero_estimate(
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  // const int nElem = alelem_ptr->get_nlocalele();
  const int nElem = locelem->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  // lassem_ptr->Assem_Estimate();
  locassem->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      const int loc_index  = locien->get_LIEN(e, i);

      for(int m=0; m<dof_mat; ++m)
        // row_index[dof_mat * i + m] = dof_mat * nbc_part->get_LID( m, loc_index ) + m;
        row_index[dof_mat * i + m] = dof_mat * nbc->get_LID( m, loc_index ) + m;
    }
    
    // MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
    //     lassem_ptr->Tangent, ADD_VALUES);
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        locassem->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  // Create a temporary zero solution vector to feed Natbc_Resis_KG
  PDNSolution * temp = new PDNSolution_NS( pnode.get(), 0, false );

  // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  NatBC_Resis_KG( 0.0, 0.1, temp, temp, ebc_part, gbc );

  delete temp;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  // for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );
  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

// void PGAssem_NS_FEM::Assem_mass_residual(
//     const PDNSolution * const &sol_a,
//     const ALocal_Elem * const &alelem_ptr,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &elementv,
//     FEAElement * const &elements,
//     FEAElement * const &elementvs,
//     const IQuadPts * const &quad_v,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const FEANode * const &fnode_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const ALocal_WeakBC * const &wbc_part )
void PGAssem_NS_FEM::Assem_mass_residual(
    const PDNSolution * const &sol_a)
{
  // const int nElem = alelem_ptr->get_nlocalele();
  const int nElem = locelem->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  double * array_a = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );

  for(int ee=0; ee<nElem; ++ee)
  {
    // lien_ptr->get_LIEN(ee, IEN_e);
    locien->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    // fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // lassem_ptr->Assem_Mass_Residual( local_a, elementv,
    //     ectrl_x, ectrl_y, ectrl_z, quad_v );
    locassem->Assem_Mass_Residual( local_a, ectrl_x, ectrl_y, ectrl_z );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        // row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
        row_index[dof_mat*ii+mm] = dof_mat * nbc -> get_LID(mm, IEN_e[ii]) + mm;
    }
    
    // MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
    //     lassem_ptr->Tangent, ADD_VALUES);
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        locassem->Tangent, ADD_VALUES);

    // VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  // Weak_EssBC_G(0, 0, sol_a, lassem_ptr, elementvs, quad_s,
  //   lien_ptr, fnode_ptr, nbc_part, wbc_part);
  Weak_EssBC_G(0, 0, sol_a);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

// void PGAssem_NS_FEM::Assem_residual(
//     const PDNSolution * const &sol_a,
//     const PDNSolution * const &sol_b,
//     const PDNSolution * const &dot_sol_np1,
//     const PDNSolution * const &sol_np1,
//     const double &curr_time,
//     const double &dt,
//     const ALocal_Elem * const &alelem_ptr,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &elementv,
//     FEAElement * const &elements,
//     FEAElement * const &elementvs,
//     const IQuadPts * const &quad_v,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const FEANode * const &fnode_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const IGenBC * const &gbc,
//     const ALocal_WeakBC * const &wbc_part )
void PGAssem_NS_FEM::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const int nElem = locelem->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for( int ee=0; ee<nElem; ++ee )
  {
    locien->get_LIEN(ee, IEN_e);
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
    
    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
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
  
  // Backflow stabilization residual contribution
  BackFlow_G( sol_b, ebc_part );

  // Resistance type boundary condition
  NatBC_Resis_G( curr_time, dt, dot_sol_np1, sol_np1, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_G( curr_time, dt, sol_b );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_G( ii );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

// void PGAssem_NS_FEM::Assem_tangent_residual(
//     const PDNSolution * const &sol_a,
//     const PDNSolution * const &sol_b,
//     const PDNSolution * const &dot_sol_np1,
//     const PDNSolution * const &sol_np1,
//     const double &curr_time,
//     const double &dt,
//     const ALocal_Elem * const &alelem_ptr,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &elementv,
//     FEAElement * const &elements,
//     FEAElement * const &elementvs,
//     const IQuadPts * const &quad_v,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const FEANode * const &fnode_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const IGenBC * const &gbc,
//     const ALocal_WeakBC * const &wbc_part )
void PGAssem_NS_FEM::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const int nElem = locelem->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

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
        row_index[dof_mat*ii + mm] = dof_mat*nbc->get_LID(mm, IEN_e[ii])+mm;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        locassem->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
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

  // Backflow stabilization residual & tangent contribution
  BackFlow_KG( dt, sol_b, ebc_part );

  // Resistance type boundary condition
  NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_KG( curr_time, dt, sol_b );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

// void PGAssem_NS_FEM::NatBC_G( const double &curr_time, const double &dt,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part )
void PGAssem_NS_FEM::NatBC_G( const double &curr_time, const double &dt,
    const ALocal_EBC * const &ebc_part )
{
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      locassem -> Assem_Residual_EBC( ebc_id, curr_time, dt,
          sctrl_x, sctrl_y, sctrl_z );

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

// void PGAssem_NS_FEM::BackFlow_G( 
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part )
void PGAssem_NS_FEM::BackFlow_G( 
  const PDNSolution * const &sol,
  const ALocal_EBC * const &ebc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  sol->GetLocalArray( array );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array, LSIEN, snLocBas, local);

      locassem->Assem_Residual_BackFlowStab( local, 
          sctrl_x, sctrl_y, sctrl_z );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc -> get_LID(mm, LSIEN[ii]) + mm;
      }

      VecSetValues(G, dof_mat*snLocBas, srow_index, locassem->sur_Residual, ADD_VALUES);
    }
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

// void PGAssem_NS_FEM::BackFlow_KG( const double &dt,
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part )
void PGAssem_NS_FEM::BackFlow_KG( const double &dt,
    const PDNSolution * const &sol,
    const ALocal_EBC * const &ebc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  sol->GetLocalArray( array );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array, LSIEN, snLocBas, local);

      locassem->Assem_Tangent_Residual_BackFlowStab( dt, local,
          sctrl_x, sctrl_y, sctrl_z );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc -> get_LID(mm, LSIEN[ii]) + mm;
      }

      MatSetValues(K, dof_mat*snLocBas, srow_index, dof_mat*snLocBas, srow_index,
          locassem->sur_Tangent, ADD_VALUES);

      VecSetValues(G, dof_mat*snLocBas, srow_index, locassem->sur_Residual, ADD_VALUES);
    }
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

// double PGAssem_NS_FEM::Assem_surface_flowrate(
//     const PDNSolution * const &vec,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_EBC * const &ebc_part,
//     const int &ebc_id )
double PGAssem_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += locassem -> get_flowrate( local, sctrl_x,
        sctrl_y, sctrl_z );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}

// double PGAssem_NS_FEM::Assem_surface_flowrate(
//     const PDNSolution * const &vec,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_InflowBC * const &infbc_part,
//     const int &nbc_id )
double PGAssem_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    const ALocal_InflowBC * const &infbc_part,
    const int &infnbc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell(infnbc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( infnbc_id, ee, LSIEN );

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( infnbc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += locassem -> get_flowrate( local, sctrl_x,
        sctrl_y, sctrl_z );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}

// double PGAssem_NS_FEM::Assem_surface_ave_pressure(
//     const PDNSolution * const &vec,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_EBC * const &ebc_part,
//     const int &ebc_id )
double PGAssem_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN );

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    locassem-> get_pressure_area( local, sctrl_x, sctrl_y,
        sctrl_z, ele_pres, ele_area );

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}

// double PGAssem_NS_FEM::Assem_surface_ave_pressure(
//     const PDNSolution * const &vec,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_InflowBC * const &infbc_part,
//     const int &nbc_id )
double PGAssem_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    const ALocal_InflowBC * const &infbc_part,
    const int &infnbc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell(infnbc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( infnbc_id, ee, LSIEN );

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( infnbc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    locassem-> get_pressure_area( local, sctrl_x, sctrl_y,
        sctrl_z, ele_pres, ele_area );

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}

// void PGAssem_NS_FEM::NatBC_Resis_G(
//     const double &curr_time, const double &dt,
//     const PDNSolution * const &dot_sol,
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const IGenBC * const &gbc )
void PGAssem_NS_FEM::NatBC_Resis_G(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id from solution vector dot_sol
    const double dot_flrate = Assem_outlet_flowrate( dot_sol, ebc_id ); 

    // Calculate flow rate for face with ebc_id from solution vector sol
    const double flrate = Assem_outlet_flowrate( sol, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate, curr_time + dt );

    // P_n+alpha_f
    // lassem_ptr->get_model_para_1() gives alpha_f 
    const double val = P_n + locassem->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // Here, val is Pressure, and is used as the surface traction h = P I 
      // to calculate the boundary integral
      locassem->Assem_Residual_EBC_Resistance(ebc_id, val,
          sctrl_x, sctrl_y, sctrl_z);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = locassem->sur_Residual[4*ii+1];
        Res[3*ii+1] = locassem->sur_Residual[4*ii+2];
        Res[3*ii+2] = locassem->sur_Residual[4*ii+3];

        srow_idx[3*ii+0] = dof_mat * nbc->get_LID(1, LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc->get_LID(2, LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc->get_LID(3, LSIEN[ii]) + 3;
      }

      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
}

// void PGAssem_NS_FEM::NatBC_Resis_KG(
//     const double &curr_time, const double &dt,
//     const PDNSolution * const &dot_sol,
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_s,
//     const IQuadPts * const &quad_s,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_EBC * const &ebc_part,
//     const IGenBC * const &gbc )
void PGAssem_NS_FEM::NatBC_Resis_KG(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const double a_f = locassem -> get_model_para_1();

  // dd_dv = dt x alpha_f x gamma
  const double dd_dv = dt * a_f * locassem->get_model_para_2();

  // Allocate the vector to hold the residual on each surface element
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];
  PetscScalar * Tan;
  PetscInt * scol_idx;
  Vector_3 out_n;
  std::vector<double> intNB;
  std::vector<int> map_Bj;

  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id and MPI_Allreduce them
    // Here, dot_sol is the solution at time step n+1 (not n+alpha_f!)
    const double dot_flrate = Assem_outlet_flowrate( dot_sol, ebc_id ); 

    // Calculate flow rate for face with ebc_id and MPI_Allreduce them
    // Here, sol is the solution at time step n+1 (not n+alpha_f!)
    const double flrate = Assem_outlet_flowrate( sol, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate, curr_time + dt );

    // P_n+alpha_f 
    const double resis_val = P_n + a_f * (P_np1 - P_n);

    // Get the (potentially approximated) m := dP/dQ
    const double m_val = gbc -> get_m( ebc_id, dot_flrate, flrate );

    // Get the (potentially approximated) n := dP/d(dot_Q)
    const double n_val = gbc -> get_n( ebc_id, dot_flrate, flrate );
    
    // Define alpha_f * n + alpha_f * gamma * dt * m
    // coef a^t a enters as the consistent tangent for the resistance-type bc
    const double coef = a_f * n_val + dd_dv * m_val;

    const int num_face_nodes = ebc_part -> get_num_face_nodes(ebc_id);
    if(num_face_nodes > 0)
    {
      Tan = new PetscScalar [snLocBas * 3 * num_face_nodes * 3];
      scol_idx = new PetscInt [num_face_nodes * 3];
      out_n  = ebc_part -> get_outvec( ebc_id );
      intNB  = ebc_part -> get_intNA( ebc_id );
      map_Bj = ebc_part -> get_LID( ebc_id );
    }
    else
    {
      Tan = nullptr;
      scol_idx = nullptr;
    }

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // For here, we scale the int_NA nx/y/z by factor 1
      locassem->Assem_Residual_EBC_Resistance(ebc_id, 1.0,
          sctrl_x, sctrl_y, sctrl_z);

      // Residual vector is scaled by the resistance value
      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = resis_val * locassem->sur_Residual[4*ii+1];
        Res[3*ii+1] = resis_val * locassem->sur_Residual[4*ii+2];
        Res[3*ii+2] = resis_val * locassem->sur_Residual[4*ii+3];
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            // Residual[4*A+ii+1] is intNB[A]*out_n[ii]
            Tan[temp_row + 3*B + 0] = coef * locassem->sur_Residual[4*A+ii+1] * intNB[B] * out_n.x();
            Tan[temp_row + 3*B + 1] = coef * locassem->sur_Residual[4*A+ii+1] * intNB[B] * out_n.y();
            Tan[temp_row + 3*B + 2] = coef * locassem->sur_Residual[4*A+ii+1] * intNB[B] * out_n.z();
          }
        }
      }

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_idx[3*ii+0] = dof_mat * nbc->get_LID(1,LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc->get_LID(2,LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc->get_LID(3,LSIEN[ii]) + 3;
      }

      for(int ii=0; ii<num_face_nodes; ++ii)
      {
        scol_idx[ii*3+0] = dof_mat * map_Bj[ii*3+0] + 1;
        scol_idx[ii*3+1] = dof_mat * map_Bj[ii*3+1] + 2;
        scol_idx[ii*3+2] = dof_mat * map_Bj[ii*3+2] + 3;
      }

      MatSetValues(K, snLocBas*3, srow_idx, num_face_nodes*3, scol_idx, Tan, ADD_VALUES);
      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }

    if( num_face_nodes > 0 ) 
    {
      delete [] Tan; Tan = nullptr;
      delete [] scol_idx; scol_idx = nullptr;
    }
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
}

// void PGAssem_NS_FEM::Weak_EssBC_KG(
//     const double &curr_time, const double &dt,
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_vs,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const FEANode * const &fnode_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_WeakBC * const &wbc_part)
void PGAssem_NS_FEM::Weak_EssBC_KG(
    const double &curr_time, const double &dt,
    const PDNSolution * const &sol )
{
  const int loc_dof {dof_mat * nLocBas};
  double * array_b = new double [nlgn * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );

  const int num_wele {wbc->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc->get_part_vol_ele_id(ee)};

    locien->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);

    fnode->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc->get_ele_face_id(ee)};

    locassem->Assem_Tangent_Residual_Weak(curr_time, dt, local_b,
      ctrl_x, ctrl_y, ctrl_z, face_id);

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc->get_LID(mm, IEN_v[ii]) + mm;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index, locassem->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

// void PGAssem_NS_FEM::Weak_EssBC_G(
//     const double &curr_time, const double &dt,
//     const PDNSolution * const &sol,
//     IPLocAssem * const &lassem_ptr,
//     FEAElement * const &element_vs,
//     const IQuadPts * const &quad_s,
//     const ALocal_IEN * const &lien_ptr,
//     const FEANode * const &fnode_ptr,
//     const ALocal_NBC * const &nbc_part,
//     const ALocal_WeakBC * const &wbc_part)
void PGAssem_NS_FEM::Weak_EssBC_G(
    const double &curr_time, const double &dt,
    const PDNSolution * const &sol )
{
  const int loc_dof {dof_mat * nLocBas};
  double * array_b = new double [nlgn * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );

  const int num_wele {wbc->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc->get_part_vol_ele_id(ee)};

    locien->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);

    fnode->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc->get_ele_face_id(ee)};
    
    locassem->Assem_Residual_Weak(curr_time, dt, local_b,
      ctrl_x, ctrl_y, ctrl_z, face_id);

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc->get_LID(mm, IEN_v[ii]) + mm;
    }

    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

// EOF
