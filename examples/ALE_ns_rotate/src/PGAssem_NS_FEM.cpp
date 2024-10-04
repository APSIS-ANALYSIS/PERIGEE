#include "PGAssem_NS_FEM.hpp"

PGAssem_NS_FEM::PGAssem_NS_FEM(
    IPLocAssem * const &locassem_ptr,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    FEAElement * const &elementvs_rotated,
    const IQuadPts * const &quads,
    IQuadPts * const &free_quad,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NBC * const &part_nbc,
    const ALocal_EBC * const &part_ebc,
    const ALocal_Interface * const &part_itf,
    SI_T::SI_solution * const &SI_sol,
    SI_T::SI_quad_point * const &SI_qp,
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat( locassem_ptr->get_dof_mat() ),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn( pnode_ptr->get_nlocghonode() ),
  snLocBas( 0 ),
  anci( locassem_ptr, elementvs, elementvs_rotated, elements,
    quads, free_quad, part_itf, SI_sol, SI_qp )
{
  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),
      "PGAssem_NS_FEM::dof_sol != locassem_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != part_nbc->get_dof_LID(),
      "PGAssem_NS_FEM::dof_mat != part_nbc->get_dof_LID(). \n");

  // Make sure that the surface element's number of local basis are 
  // the same. This is an assumption in this assembly routine.
  if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_NS_FEM, snLocBas has to be uniform. \n");
  }

  const int nlocrow = dof_mat * pnode_ptr->get_nlocalnode();

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

  SI_qp->search_all_opposite_point(elementvs, elementvs_rotated, elements, quads, free_quad, part_itf, SI_sol);

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, 
      elements, elementvs, elementvs_rotated, quads, free_quad, aien_ptr, pnode_ptr, part_nbc, part_ebc, part_itf, SI_sol, SI_qp, gbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation
 
  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

  // Assem_L2_proj_mat(elementvs, elements, quads, part_itf, SI_sol);
}

PGAssem_NS_FEM::~PGAssem_NS_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
  VecDestroy(&L2_proj_lhs);
  VecDestroy(&L2_proj_sol);
  VecDestroy(&L2_proj_sol_x);
  VecDestroy(&L2_proj_sol_y);
  VecDestroy(&L2_proj_sol_z);
  VecDestroy(&L2_proj_mvelo);
  MatDestroy(&L2_proj_mat);
}

void PGAssem_NS_FEM::EssBC_KG(
    const ALocal_NBC * const &nbc_part, const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc_part->get_LDN(field, i) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(field);
  if(local_sla > 0)
  {
    for(int i=0; i<local_sla; ++i)
    {
      const int row = nbc_part->get_LPSN(field, i) * dof_mat + field;
      const int col = nbc_part->get_LPMN(field, i) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_NS_FEM::EssBC_G( const ALocal_NBC * const &nbc_part, 
    const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc_part->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc_part->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_NS_FEM::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_Interface * const &itf_part,
        const SI_T::SI_solution * const &SI_sol,
        const SI_T::SI_quad_point * const &SI_qp,
        const IGenBC * const &gbc )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  lassem_ptr->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      const int loc_index  = lien_ptr->get_LIEN(e, i);

      for(int m=0; m<dof_mat; ++m)
        row_index[dof_mat * i + m] = dof_mat * nbc_part->get_LID( m, loc_index ) + m;
    }
    
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  // Create a temporary zero solution vector to feed Natbc_Resis_KG
  PDNSolution * temp = new PDNSolution_NS( node_ptr, 0, false );

  // // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  // NatBC_Resis_KG( 0.0, 0.1, temp, temp, lassem_ptr, elements, quad_s, nbc_part, ebc_part, gbc );

  delete temp;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_NS_FEM::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &mdisp,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    FEAElement * const &elementvs_rotated,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IQuadPts * const &free_quad,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_WeakBC * const &wbc_part,
    const ALocal_Interface * const &itf_part,
    const SI_T::SI_solution * const &SI_sol,
    const SI_T::SI_quad_point * const &SI_qp )
{
  const int nElem = alelem_ptr->get_nlocalele();
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
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Mass_Residual( local_a, elementv,
        ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
    }
    
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
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
  // Default mdisp = mvelo = zeros at t = 0.
  Weak_EssBC_G(0, 0, sol_a, mdisp, mdisp, alelem_ptr, lassem_ptr, elementvs, quad_s,
    lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // Surface integral from Nitsche method
  Interface_G(0, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part, SI_sol, SI_qp);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_NS_FEM::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &mvelo,
    const PDNSolution * const &mdisp,    
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    FEAElement * const &elementvs_rotated,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IQuadPts * const &free_quad,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const ALocal_Interface * const &itf_part,
    const SI_T::SI_solution * const &SI_sol,
    const SI_T::SI_quad_point * const &SI_qp )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * array_mvelo = new double [nlgn * 3];
  double * array_mdisp = new double [nlgn * 3];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  double * local_mvelo = new double [nLocBas * 3];
  double * local_mdisp = new double [nLocBas * 3];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );
  mvelo->GetLocalArray( array_mvelo );
  mdisp->GetLocalArray( array_mdisp );

  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);
    GetLocal(array_mvelo, 3, IEN_e, local_mvelo);  // 3 is the DOFs of mvelo/mdisp
    GetLocal(array_mdisp, 3, IEN_e, local_mdisp);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b, local_mvelo, local_mdisp,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
    }

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] array_mvelo; array_mvelo = nullptr;
  delete [] array_mdisp; array_mdisp = nullptr;
  delete [] local_mvelo; local_mvelo = nullptr;
  delete [] local_mdisp; local_mdisp = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
   
  // Backflow stabilization residual contribution
  // BackFlow_G( sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_G( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
  //    nbc_part, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_G(curr_time, dt, sol_b, mvelo, mdisp, alelem_ptr, lassem_ptr, elementvs, quad_s,
      lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // // For Poiseuille flow
  // NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // Surface integral from Nitsche method
  Interface_G(dt, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part, SI_sol, SI_qp);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_G( nbc_part, ii );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_NS_FEM::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &mvelo,
    const PDNSolution * const &mdisp, 
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    FEAElement * const &elementvs_rotated,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IQuadPts * const &free_quad,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const ALocal_Interface * const &itf_part,
    const SI_T::SI_solution * const &SI_sol,
    const SI_T::SI_quad_point * const &SI_qp )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * array_mvelo = new double [nlgn * 3];
  double * array_mdisp = new double [nlgn * 3];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  double * local_mvelo = new double [nLocBas * 3];
  double * local_mdisp = new double [nLocBas * 3];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );
  mvelo->GetLocalArray( array_mvelo );
  mdisp->GetLocalArray( array_mdisp );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);
    GetLocal(array_mvelo, 3, IEN_e, local_mvelo);  // 3 is the DOFs of mvelo/mdisp
    GetLocal(array_mdisp, 3, IEN_e, local_mdisp);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b, local_mvelo, local_mdisp,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_e[ii])+mm;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
      lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] array_mvelo; array_mvelo = nullptr;
  delete [] array_mdisp; array_mdisp = nullptr;
  delete [] local_mvelo; local_mvelo = nullptr;
  delete [] local_mdisp; local_mdisp = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Backflow stabilization residual & tangent contribution
  // BackFlow_KG( dt, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
  //    nbc_part, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_KG(curr_time, dt, sol_b, mvelo, mdisp, alelem_ptr, lassem_ptr, elementvs, quad_s,
    lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // // For Poiseuille flow
  // NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // Surface integral from Nitsche method (Only assemble G at present)
  Interface_G(dt, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part, SI_sol, SI_qp);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_NS_FEM::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
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

      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
    }
  }

  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

void PGAssem_NS_FEM::BackFlow_G( 
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [dof_sol * snLocBas];
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

      lassem_ptr->Assem_Residual_BackFlowStab( local, element_s,
        sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
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

void PGAssem_NS_FEM::BackFlow_KG( const double &dt,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [dof_sol * snLocBas];
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

      lassem_ptr->Assem_Tangent_Residual_BackFlowStab( dt, local,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      MatSetValues(K, dof_mat*snLocBas, srow_index, dof_mat*snLocBas, srow_index,
          lassem_ptr->sur_Tangent, ADD_VALUES);

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
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

double PGAssem_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
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

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
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

double PGAssem_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_InflowBC * const &infbc_part,
    const int &nbc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell(nbc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( nbc_id, ee, LSIEN );

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( nbc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
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

double PGAssem_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
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
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

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

double PGAssem_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_InflowBC * const &infbc_part,
    const int &nbc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell(nbc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( nbc_id, ee, LSIEN );

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( nbc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

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

void PGAssem_NS_FEM::NatBC_Resis_G(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
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
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id from solution vector sol
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate, curr_time + dt );

    // P_n+alpha_f
    // lassem_ptr->get_model_para_1() gives alpha_f 
    const double val = P_n + lassem_ptr->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // Here, val is Pressure, and is used as the surface traction h = P I 
      // to calculate the boundary integral
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, val,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = lassem_ptr->Residual[4*ii+1];
        Res[3*ii+1] = lassem_ptr->Residual[4*ii+2];
        Res[3*ii+2] = lassem_ptr->Residual[4*ii+3];

        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1, LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2, LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3, LSIEN[ii]) + 3;
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

void PGAssem_NS_FEM::NatBC_Resis_KG(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const double a_f = lassem_ptr -> get_model_para_1();

  // dd_dv = dt x alpha_f x gamma
  const double dd_dv = dt * a_f * lassem_ptr->get_model_para_2();

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
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id and MPI_Allreduce them
    // Here, sol is the solution at time step n+1 (not n+alpha_f!)
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

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
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, 1.0,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      // Residual vector is scaled by the resistance value
      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = resis_val * lassem_ptr->Residual[4*ii+1];
        Res[3*ii+1] = resis_val * lassem_ptr->Residual[4*ii+2];
        Res[3*ii+2] = resis_val * lassem_ptr->Residual[4*ii+3];
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            // Residual[4*A+ii+1] is intNB[A]*out_n[ii]
            Tan[temp_row + 3*B + 0] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.x();
            Tan[temp_row + 3*B + 1] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.y();
            Tan[temp_row + 3*B + 2] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.z();
          }
        }
      }

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1,LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2,LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3,LSIEN[ii]) + 3;
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

void PGAssem_NS_FEM::Weak_EssBC_KG(
    const double &curr_time, const double &dt,
    const PDNSolution * const &sol,
    const PDNSolution * const &mvelo,
    const PDNSolution * const &mdisp,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_vs,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_WeakBC * const &wbc_part)
{
  const int loc_dof {dof_mat * nLocBas};
  double * array_b = new double [nlgn * dof_sol];
  double * array_mvelo = new double [nlgn * 3];
  double * array_mdisp = new double [nlgn * 3];
  double * local_b = new double [nLocBas * dof_sol];
  double * local_mvelo = new double [nLocBas * 3];
  double * local_mdisp = new double [nLocBas * 3];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );
  mvelo->GetLocalArray( array_mvelo );
  mdisp->GetLocalArray( array_mdisp );

  const int num_wele {wbc_part->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc_part->get_part_vol_ele_id(ee)};

    lien_ptr->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);
    GetLocal(array_mvelo, 3, IEN_v, local_mvelo);
    GetLocal(array_mdisp, 3, IEN_v, local_mdisp);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc_part->get_ele_face_id(ee)};

    lassem_ptr->Assem_Tangent_Residual_Weak(curr_time, dt, local_b, local_mvelo, local_mdisp,
      element_vs, ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_v[ii]) + mm;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index, lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] array_mvelo; array_mvelo = nullptr;
  delete [] array_mdisp; array_mdisp = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] local_mvelo; local_mvelo = nullptr;
  delete [] local_mdisp; local_mdisp = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

void PGAssem_NS_FEM::Weak_EssBC_G(
    const double &curr_time, const double &dt,
    const PDNSolution * const &sol,
    const PDNSolution * const &mvelo,
    const PDNSolution * const &mdisp,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_vs,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_WeakBC * const &wbc_part)
{
  const int loc_dof {dof_mat * nLocBas};
  double * array_b = new double [nlgn * dof_sol];
  double * array_mvelo = new double [nlgn * 3];
  double * array_mdisp = new double [nlgn * 3];
  double * local_b = new double [nLocBas * dof_sol];
  double * local_mvelo = new double [nLocBas * 3];
  double * local_mdisp = new double [nLocBas * 3];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );
  mvelo->GetLocalArray( array_mvelo );
  mdisp->GetLocalArray( array_mdisp );

  const int num_wele {wbc_part->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc_part->get_part_vol_ele_id(ee)};

    lien_ptr->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);
    GetLocal(array_mvelo, 3, IEN_v, local_mvelo);
    GetLocal(array_mdisp, 3, IEN_v, local_mdisp);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc_part->get_ele_face_id(ee)};
    
    lassem_ptr->Assem_Residual_Weak(curr_time, dt, local_b, local_mvelo, local_mdisp,
      element_vs, ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_v[ii]) + mm;
    }

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] array_mvelo; array_mvelo = nullptr;
  delete [] array_mdisp; array_mdisp = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] local_mvelo; local_mvelo = nullptr;
  delete [] local_mdisp; local_mdisp = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

void PGAssem_NS_FEM::Interface_G(
  const double &dt,
  IPLocAssem * const &lassem_ptr,
  FEAElement * const &fixed_elementv,
  FEAElement * const &rotated_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  IQuadPts * const &free_quad,
  const ALocal_Interface * const &itf_part,
  const SI_T::SI_solution * const &SI_sol,
  const SI_T::SI_quad_point * const &SI_qp )
{
  anci.A_dt = dt;
  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  double * curPt_x = new double [nLocBas];
  double * curPt_y = new double [nLocBas];
  double * curPt_z = new double [nLocBas];

  int * fixed_local_ien = new int [nLocBas];
  double * fixed_local_sol = new double [nLocBas * dof_sol];

  int * rotated_local_ien = new int [nLocBas];
  double * proj_rotated_local_sol = new double [nLocBas * dof_sol];
  double * proj_rotated_local_sol_x = new double [nLocBas * 3];
  double * proj_rotated_local_sol_y = new double [nLocBas * 3];
  double * proj_rotated_local_sol_z = new double [nLocBas * 3];
  double * proj_rotated_local_mvelo = new double [nLocBas * 3];
  double * rotated_local_disp = new double [nLocBas * 3];

  PetscInt * fixed_row_index = new PetscInt [nLocBas * dof_mat];
  PetscInt * rotated_row_index = new PetscInt [nLocBas * dof_mat];

  const int num_itf {itf_part->get_num_itf()};

  int ele_tag {-1};
  int rotated_ee {-1};
  double rotated_xi {1.0 / 3.0};
  double rotated_eta {1.0 / 3.0};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    // SYS_T::commPrint("itf_id = %d\n", itf_id);
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    const int fixed_face_nqp {quad_s->get_num_quadPts()};

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      // SYS_T::commPrint("  fixed_ee = %d\n", ee);
      // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      std::vector<double> R(nLocBas, 0.0);

      // Get the local ien and local sol of this fixed element
      SI_sol->get_fixed_local(itf_part, itf_id, ee, fixed_local_ien, fixed_local_sol);

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        fixed_elementv->get_R(qua, &R[0]);

        SI_qp->get_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi, rotated_eta);

        const int rotated_face_id {itf_part->get_rotated_face_id(itf_id, ele_tag, rotated_ee)};

        itf_part->get_rotated_ele_ctrlPts(itf_id, ele_tag, rotated_ee, ctrl_x, ctrl_y, ctrl_z);

        SI_sol->get_rotated_mdisp(itf_part, itf_id, ele_tag, rotated_ee, rotated_local_ien, rotated_local_disp);
        get_currPts(ctrl_x, ctrl_y, ctrl_z, rotated_local_disp, nLocBas, curPt_x, curPt_y, curPt_z);

        free_quad->set_qp( rotated_xi, rotated_eta );

        rotated_elementv->buildBasis(rotated_face_id, free_quad, curPt_x, curPt_y, curPt_z);

        SI_sol->get_projected_rotated_local(itf_id, fixed_local_ien, proj_rotated_local_sol,
          proj_rotated_local_sol_x, proj_rotated_local_sol_y, proj_rotated_local_sol_z, proj_rotated_local_mvelo);

        const double qw = quad_s->get_qw(qua);

        lassem_ptr->Assem_Residual_itf(qua, qw, dt, fixed_elementv, rotated_elementv, fixed_local_sol,
          proj_rotated_local_sol, proj_rotated_local_sol_x, proj_rotated_local_sol_y, proj_rotated_local_sol_z,
          proj_rotated_local_mvelo);

        for(int ii{0}; ii < nLocBas; ++ii)
        {
          for(int mm{0}; mm < dof_mat; ++mm)
          {
            fixed_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_fixed_LID(itf_id, mm, fixed_local_ien[ii]) + mm;
            rotated_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_rotated_LID(itf_id, mm, rotated_local_ien[ii]) + mm;
          }
        }

        VecSetValues(G, loc_dof, fixed_row_index, lassem_ptr->Residual_s, ADD_VALUES);
        VecSetValues(G, loc_dof, rotated_row_index, lassem_ptr->Residual_r, ADD_VALUES);
      }
    }
  }

  delete [] fixed_local_ien; fixed_local_ien = nullptr;
  delete [] fixed_local_sol; fixed_local_sol = nullptr;
  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] proj_rotated_local_sol; proj_rotated_local_sol = nullptr;
  delete [] proj_rotated_local_sol_x; proj_rotated_local_sol_x = nullptr;
  delete [] proj_rotated_local_sol_y; proj_rotated_local_sol_y = nullptr;
  delete [] proj_rotated_local_sol_z; proj_rotated_local_sol_z = nullptr;
  delete [] proj_rotated_local_mvelo; proj_rotated_local_mvelo = nullptr;
  delete [] rotated_local_disp; rotated_local_disp = nullptr;

  delete [] fixed_row_index; fixed_row_index = nullptr;
  delete [] rotated_row_index; rotated_row_index = nullptr;

  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;

  delete [] curPt_x; curPt_x = nullptr;
  delete [] curPt_y; curPt_y = nullptr;
  delete [] curPt_z; curPt_z = nullptr;
}

void PGAssem_NS_FEM::Interface_K_MF(Vec &X, Vec &Y)
{
  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  double * curPt_x = new double [nLocBas];
  double * curPt_y = new double [nLocBas];
  double * curPt_z = new double [nLocBas];

  int * fixed_local_ien = new int [nLocBas];
  double * fixed_local_sol = new double [nLocBas * dof_sol];

  int * rotated_local_ien = new int [nLocBas];
  double * proj_rotated_local_sol = new double [nLocBas * dof_sol];
  double * proj_rotated_local_sol_x = new double [nLocBas * 3];
  double * proj_rotated_local_sol_y = new double [nLocBas * 3];
  double * proj_rotated_local_sol_z = new double [nLocBas * 3];
  double * proj_rotated_local_mvelo = new double [nLocBas * 3];
  double * rotated_local_disp = new double [nLocBas * 3];

  PetscInt * fixed_row_index = new PetscInt [nLocBas * dof_mat];
  PetscInt * rotated_row_index = new PetscInt [nLocBas * dof_mat];

  const int num_itf {anci.A_itf_part->get_num_itf()};

  int ele_tag {-1};
  int rotated_ee {-1};
  double rotated_xi {1.0 / 3.0};
  double rotated_eta {1.0 / 3.0};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    // SYS_T::commPrint("itf_id = %d\n", itf_id);
    const int num_fixed_elem = anci.A_itf_part->get_num_fixed_ele(itf_id);

    const int fixed_face_nqp {anci.A_quad_s->get_num_quadPts()};

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      // SYS_T::commPrint("  fixed_ee = %d\n", ee);
      // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

      anci.A_itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {anci.A_itf_part->get_fixed_face_id(itf_id, ee)};

      anci.A_fixed_elementv->buildBasis(fixed_face_id, anci.A_quad_s, ctrl_x, ctrl_y, ctrl_z);

      std::vector<double> R(nLocBas, 0.0);

      // Get the local ien and local sol of this fixed element
      anci.A_SI_sol->get_fixed_local(anci.A_itf_part, itf_id, ee, fixed_local_ien, fixed_local_sol);

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        anci.A_fixed_elementv->get_R(qua, &R[0]);

        anci.A_SI_qp->get_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi, rotated_eta);

        const int rotated_face_id {anci.A_itf_part->get_rotated_face_id(itf_id, ele_tag, rotated_ee)};

        anci.A_itf_part->get_rotated_ele_ctrlPts(itf_id, ele_tag, rotated_ee, ctrl_x, ctrl_y, ctrl_z);

        anci.A_SI_sol->get_rotated_mdisp(anci.A_itf_part, itf_id, ele_tag, rotated_ee, rotated_local_ien, rotated_local_disp);
        get_currPts(ctrl_x, ctrl_y, ctrl_z, rotated_local_disp, nLocBas, curPt_x, curPt_y, curPt_z);

        anci.A_free_quad->set_qp( rotated_xi, rotated_eta );

        anci.A_rotated_elementv->buildBasis(rotated_face_id, anci.A_free_quad, curPt_x, curPt_y, curPt_z);

        anci.A_SI_sol->get_projected_rotated_local(itf_id, fixed_local_ien, proj_rotated_local_sol,
          proj_rotated_local_sol_x, proj_rotated_local_sol_y, proj_rotated_local_sol_z, proj_rotated_local_mvelo);

        const double qw = anci.A_quad_s->get_qw(qua);

        anci.A_lassemptr->Assem_Tangent_itf_MF(qua, qw, anci.A_dt, anci.A_fixed_elementv, anci.A_rotated_elementv,
          fixed_local_sol, proj_rotated_local_sol, proj_rotated_local_mvelo);

        for(int ii{0}; ii < nLocBas; ++ii)
        {
          for(int mm{0}; mm < dof_mat; ++mm)
          {
            fixed_row_index[dof_mat * ii + mm] = dof_mat * anci.A_itf_part->get_fixed_LID(itf_id, mm, fixed_local_ien[ii]) + mm;
            rotated_row_index[dof_mat * ii + mm] = dof_mat * anci.A_itf_part->get_rotated_LID(itf_id, mm, rotated_local_ien[ii]) + mm;
          }
        }

        // Matrix-free method
        local_MatMult_MF(loc_dof, fixed_row_index, fixed_row_index, anci.A_lassemptr->Tangent_ss, X, Y);
        local_MatMult_MF(loc_dof, rotated_row_index, rotated_row_index, anci.A_lassemptr->Tangent_rr, X, Y);
        local_MatMult_MF(loc_dof, fixed_row_index, rotated_row_index, anci.A_lassemptr->Tangent_sr, X, Y);
        local_MatMult_MF(loc_dof, rotated_row_index, fixed_row_index, anci.A_lassemptr->Tangent_rs, X, Y);
      }
    }

    // Synchronize all CPUs, see PETSc_Tool.hpp: PETSc_T::Scatter
    const int max_num_fixed_elem = anci.A_itf_part->get_max_num_fixed_ele(itf_id);
    
    double * local_X = new double[loc_dof];
    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        fixed_row_index[dof_mat * ii + mm] = 0;
    }

    for(int ee{0}; ee < max_num_fixed_elem - num_fixed_elem; ++ee)
    {
      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        PETSc_T::Scatter(X, fixed_row_index, loc_dof, local_X);
        PETSc_T::Scatter(X, fixed_row_index, loc_dof, local_X);
        PETSc_T::Scatter(X, fixed_row_index, loc_dof, local_X);
        PETSc_T::Scatter(X, fixed_row_index, loc_dof, local_X);
      }
    }
    delete [] local_X; local_X = nullptr;
  }

  VecAssemblyBegin(Y);
  VecAssemblyEnd(Y);

  delete [] fixed_local_ien; fixed_local_ien = nullptr;
  delete [] fixed_local_sol; fixed_local_sol = nullptr;
  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] proj_rotated_local_sol; proj_rotated_local_sol = nullptr;
  delete [] proj_rotated_local_sol_x; proj_rotated_local_sol_x = nullptr;
  delete [] proj_rotated_local_sol_y; proj_rotated_local_sol_y = nullptr;
  delete [] proj_rotated_local_sol_z; proj_rotated_local_sol_z = nullptr;
  delete [] proj_rotated_local_mvelo; proj_rotated_local_mvelo = nullptr;
  delete [] rotated_local_disp; rotated_local_disp = nullptr;

  delete [] fixed_row_index; fixed_row_index = nullptr;
  delete [] rotated_row_index; rotated_row_index = nullptr;

  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;

  delete [] curPt_x; curPt_x = nullptr;
  delete [] curPt_y; curPt_y = nullptr;
  delete [] curPt_z; curPt_z = nullptr;
}

void PGAssem_NS_FEM::local_MatMult_MF(
  const int &dof,
  PetscInt * &row_index, PetscInt * &col_index,
  PetscScalar * &Mat, Vec &X, Vec &Y )
{
  double * local_X = new double [dof];
  double * local_dY = new double [dof];

  PETSc_T::Scatter(X, col_index, dof, local_X);

  for(int AA{0}; AA<dof; ++AA)
  {
    local_dY[AA] = 0.0;

    if(row_index[AA] >= 0)
    {
      for(int BB{0}; BB<dof; ++BB)
      {
        if(col_index[BB] >= 0)
          local_dY[AA] += Mat[dof * AA + BB] * local_X[BB];
      }
    }
  }

  VecSetValues(Y, dof, row_index, local_dY, ADD_VALUES);

  delete [] local_X; local_X = nullptr;
  delete [] local_dY; local_dY = nullptr;
}

void PGAssem_NS_FEM::Assem_L2_proj_mat(
  FEAElement * const &fixed_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  const ALocal_Interface * const &itf_part,
  const SI_T::SI_solution * const &SI_sol )
{
  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  PetscScalar * mat_L2_proj = new PetscScalar[loc_dof * loc_dof];
  PetscInt * row_index = new PetscInt[loc_dof];

  const int num_itf {itf_part->get_num_itf()};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    const int fixed_face_nqp {quad_s->get_num_quadPts()};

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      for(int dd=0; dd<loc_dof * loc_dof; ++dd)
        mat_L2_proj[dd] = 0.0;

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        std::vector<double> R = fixed_elementv->get_R(qua);

        const double gwts = fixed_elementv->get_detJac(qua) * quad_s->get_qw(qua);

        for(int A=0; A<nLocBas; ++A)
        {
          const double NA = R[A];
          for(int B=0; B<nLocBas; ++B)
          {
            const double NB = R[B];

            mat_L2_proj[16*nLocBas*A+4*B] += gwts * NA * NB;

            mat_L2_proj[4*nLocBas*(4*A+1)+4*B+1] += gwts * NA * NB;

            mat_L2_proj[4*nLocBas*(4*A+2)+4*B+2] += gwts * NA * NB;

            mat_L2_proj[4*nLocBas*(4*A+3)+4*B+3] += gwts * NA * NB;
          }
        }
      }

      for(int ii=0; ii<nLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
        {
          const int node = itf_part->get_fixed_lien(itf_id, nLocBas * ee + ii);
          row_index[dof_mat * ii + mm] = dof_mat * SI_sol->get_L2_proj_fixed_node_pos(itf_id, node) + mm;
        }
      }

      MatSetValues(L2_proj_mat, loc_dof, row_index, loc_dof, row_index, mat_L2_proj, ADD_VALUES);
    }
  }

  int num_inner = SI_sol->get_num_all_fixed_inner_node();
  for(int ii=0; ii<num_inner; ++ii)
  {
    for(int mm=0; mm<dof_mat; ++mm)
    {
      const int row = dof_mat * SI_sol->get_fixed_inner_node(ii) + mm;
      MatSetValue(L2_proj_mat, row, row, 1.0, ADD_VALUES);
    }
  }

  MatAssemblyBegin(L2_proj_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(L2_proj_mat, MAT_FINAL_ASSEMBLY);

  delete [] row_index; row_index = nullptr;
  delete [] mat_L2_proj; mat_L2_proj = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
}

void PGAssem_NS_FEM::Assem_L2_proj_rhs(
  FEAElement * const &fixed_elementv,
  FEAElement * const &rotated_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  IQuadPts * const &free_quad,
  const ALocal_Interface * const &itf_part,
  const SI_T::SI_solution * const &SI_sol,
  const SI_T::SI_quad_point * const &SI_qp)
{
  Clear_L2_proj_rhs();

  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  double * curPt_x = new double [nLocBas];
  double * curPt_y = new double [nLocBas];
  double * curPt_z = new double [nLocBas];

  int * rotated_local_ien = new int [nLocBas];
  double * rotated_local_sol = new double [nLocBas * dof_sol];
  double * rotated_local_mvelo = new double [nLocBas * 3];
  double * rotated_local_disp = new double [nLocBas * 3];

  PetscScalar * res_L2_proj_sol = new PetscScalar[loc_dof];
  PetscScalar * res_L2_proj_sol_x = new PetscScalar[loc_dof];
  PetscScalar * res_L2_proj_sol_y = new PetscScalar[loc_dof];
  PetscScalar * res_L2_proj_sol_z = new PetscScalar[loc_dof];
  PetscScalar * res_L2_proj_mvelo = new PetscScalar[loc_dof];
  PetscInt * row_index = new PetscInt[loc_dof];

  const int num_itf {itf_part->get_num_itf()};

  int ele_tag {-1};
  int rotated_ee {-1};
  double rotated_xi {1.0 / 3.0};
  double rotated_eta {1.0 / 3.0};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    const int fixed_face_nqp {quad_s->get_num_quadPts()};

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      for(int dd=0; dd<loc_dof; ++dd)
      {
        res_L2_proj_sol[dd] = 0.0;
        res_L2_proj_sol_x[dd] = 0.0;
        res_L2_proj_sol_y[dd] = 0.0;
        res_L2_proj_sol_z[dd] = 0.0;
        res_L2_proj_mvelo[dd] = 0.0;
      }

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        SI_qp->get_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi, rotated_eta);

        const int rotated_face_id {itf_part->get_rotated_face_id(itf_id, ele_tag, rotated_ee)};

        itf_part->get_rotated_ele_ctrlPts(itf_id, ele_tag, rotated_ee, ctrl_x, ctrl_y, ctrl_z);

        SI_sol->get_rotated_mdisp(itf_part, itf_id, ele_tag, rotated_ee, rotated_local_ien, rotated_local_disp);
        get_currPts(ctrl_x, ctrl_y, ctrl_z, rotated_local_disp, nLocBas, curPt_x, curPt_y, curPt_z);

        free_quad->set_qp( rotated_xi, rotated_eta );

        rotated_elementv->buildBasis(rotated_face_id, free_quad, curPt_x, curPt_y, curPt_z);

        SI_sol->get_rotated_local(itf_id, rotated_local_ien, rotated_local_sol, rotated_local_mvelo);

        const double gwts = fixed_elementv->get_detJac(qua) * quad_s->get_qw(qua);

        std::vector<double> Ns = fixed_elementv->get_R(qua);

        std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

        rotated_elementv -> get_R_gradR( 0, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

        double pr {0.0};
        double ur {0.0}, ur_x {0.0}, ur_y {0.0}, ur_z {0.0};
        double vr {0.0}, vr_x {0.0}, vr_y {0.0}, vr_z {0.0};
        double wr {0.0}, wr_x {0.0}, wr_y {0.0}, wr_z {0.0};
        double mur {0.0}, mvr {0.0}, mwr{0.0};

        for(int ii=0; ii<nLocBas; ++ii)
        {
          const int ii4{4 * ii};
          const int ii3{3 * ii};

          pr += rotated_local_sol[ii4 + 0] * Nr[ii];
          ur += rotated_local_sol[ii4 + 1] * Nr[ii];
          vr += rotated_local_sol[ii4 + 2] * Nr[ii];
          wr += rotated_local_sol[ii4 + 3] * Nr[ii];

          mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
          mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
          mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];

          ur_x += rotated_local_sol[ii4 + 1] * dNr_dx[ii];
          ur_y += rotated_local_sol[ii4 + 1] * dNr_dy[ii];
          ur_z += rotated_local_sol[ii4 + 1] * dNr_dz[ii];

          vr_x += rotated_local_sol[ii4 + 2] * dNr_dx[ii];
          vr_y += rotated_local_sol[ii4 + 2] * dNr_dy[ii];
          vr_z += rotated_local_sol[ii4 + 2] * dNr_dz[ii];

          wr_x += rotated_local_sol[ii4 + 3] * dNr_dx[ii];
          wr_y += rotated_local_sol[ii4 + 3] * dNr_dy[ii];
          wr_z += rotated_local_sol[ii4 + 3] * dNr_dz[ii];
        }

        for(int A=0; A<nLocBas; ++A)
        {
          const double NA = Ns[A];
          const int A4 = 4 * A;

          res_L2_proj_sol[A4] += gwts * NA * pr;
          res_L2_proj_sol[A4+1] += gwts * NA * ur;
          res_L2_proj_sol[A4+2] += gwts * NA * vr;
          res_L2_proj_sol[A4+3] += gwts * NA * wr;

          res_L2_proj_sol_x[A4+1] += gwts * NA * ur_x;
          res_L2_proj_sol_x[A4+2] += gwts * NA * vr_x;
          res_L2_proj_sol_x[A4+3] += gwts * NA * wr_x;

          res_L2_proj_sol_y[A4+1] += gwts * NA * ur_y;
          res_L2_proj_sol_y[A4+2] += gwts * NA * vr_y;
          res_L2_proj_sol_y[A4+3] += gwts * NA * wr_y;

          res_L2_proj_sol_z[A4+1] += gwts * NA * ur_z;
          res_L2_proj_sol_z[A4+2] += gwts * NA * vr_z;
          res_L2_proj_sol_z[A4+3] += gwts * NA * wr_z;

          res_L2_proj_mvelo[A4+1] += gwts * NA * mur;
          res_L2_proj_mvelo[A4+2] += gwts * NA * mvr;
          res_L2_proj_mvelo[A4+3] += gwts * NA * mwr;
        }
      }

      for(int ii=0; ii<nLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
        {
          const int node = itf_part->get_fixed_lien(itf_id, nLocBas * ee + ii);
          row_index[dof_mat * ii + mm] = dof_mat * SI_sol->get_L2_proj_fixed_node_pos(itf_id, node) + mm;
        }
      }

      VecSetValues(L2_proj_sol, loc_dof, row_index, res_L2_proj_sol, ADD_VALUES);
      VecSetValues(L2_proj_sol_x, loc_dof, row_index, res_L2_proj_sol_x, ADD_VALUES);
      VecSetValues(L2_proj_sol_y, loc_dof, row_index, res_L2_proj_sol_y, ADD_VALUES);
      VecSetValues(L2_proj_sol_z, loc_dof, row_index, res_L2_proj_sol_z, ADD_VALUES);
      VecSetValues(L2_proj_mvelo, loc_dof, row_index, res_L2_proj_mvelo, ADD_VALUES);
    }
  }

  VecAssemblyBegin(L2_proj_sol);
  VecAssemblyEnd(L2_proj_sol);

  VecAssemblyBegin(L2_proj_sol_x);
  VecAssemblyEnd(L2_proj_sol_x);

  VecAssemblyBegin(L2_proj_sol_y);
  VecAssemblyEnd(L2_proj_sol_y);

  VecAssemblyBegin(L2_proj_sol_z);
  VecAssemblyEnd(L2_proj_sol_z);

  VecAssemblyBegin(L2_proj_mvelo);
  VecAssemblyEnd(L2_proj_mvelo);

  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] rotated_local_sol; rotated_local_sol = nullptr;
  delete [] rotated_local_mvelo; rotated_local_mvelo = nullptr;
  delete [] rotated_local_disp; rotated_local_disp = nullptr;

  delete [] row_index; row_index = nullptr;
  delete [] res_L2_proj_sol; res_L2_proj_sol = nullptr;
  delete [] res_L2_proj_sol_x; res_L2_proj_sol_x = nullptr;
  delete [] res_L2_proj_sol_y; res_L2_proj_sol_y = nullptr;
  delete [] res_L2_proj_sol_z; res_L2_proj_sol_z = nullptr;
  delete [] res_L2_proj_mvelo; res_L2_proj_mvelo = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
}

void PGAssem_NS_FEM::Solve_L2_proj(PLinear_Solver_PETSc * const &lsolver_ptr)
{
  Assem_L2_proj_rhs(anci.A_fixed_elementv, anci.A_rotated_elementv, anci.A_elements,
    anci.A_quad_s, anci.A_free_quad, anci.A_itf_part, anci.A_SI_sol, anci.A_SI_qp);

  lsolver_ptr->Solve(L2_proj_mat, L2_proj_sol, L2_proj_lhs, true);
  anci.A_SI_sol->update_L2_proj_result(L2_proj_lhs, 0);


  lsolver_ptr->Solve(L2_proj_mat, L2_proj_sol_x, L2_proj_lhs, true);
  anci.A_SI_sol->update_L2_proj_result(L2_proj_lhs, 1);

  lsolver_ptr->Solve(L2_proj_mat, L2_proj_sol_y, L2_proj_lhs, true);
  anci.A_SI_sol->update_L2_proj_result(L2_proj_lhs, 2);

  lsolver_ptr->Solve(L2_proj_mat, L2_proj_sol_z, L2_proj_lhs, true);
  anci.A_SI_sol->update_L2_proj_result(L2_proj_lhs, 3);

  lsolver_ptr->Solve(L2_proj_mat, L2_proj_mvelo, L2_proj_lhs, true);
  anci.A_SI_sol->update_L2_proj_result(L2_proj_lhs, 4);

  SYS_T::commPrint("\n");
}

void PGAssem_NS_FEM::Init_L2_proj()
{
  const int n_L2proj_totalrow = 4 * anci.A_SI_sol->get_num_all_fixed_node();
  MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
    n_L2proj_totalrow, n_L2proj_totalrow,
    PETSC_DECIDE, NULL, PETSC_DECIDE, NULL, &L2_proj_mat);

  MatZeroEntries(L2_proj_mat);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_lhs);
  VecSetSizes(L2_proj_lhs, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_lhs);
  VecSet(L2_proj_lhs, 0.0);
  VecSetOption(L2_proj_lhs, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  VecCreate(PETSC_COMM_WORLD, &L2_proj_sol);
  VecSetSizes(L2_proj_sol, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_sol);
  VecSet(L2_proj_sol, 0.0);
  VecSetOption(L2_proj_sol, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_sol_x);
  VecSetSizes(L2_proj_sol_x, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_sol_x);
  VecSet(L2_proj_sol_x, 0.0);
  VecSetOption(L2_proj_sol_x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_sol_y);
  VecSetSizes(L2_proj_sol_y, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_sol_y);
  VecSet(L2_proj_sol_y, 0.0);
  VecSetOption(L2_proj_sol_y, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_sol_z);
  VecSetSizes(L2_proj_sol_z, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_sol_z);
  VecSet(L2_proj_sol_z, 0.0);
  VecSetOption(L2_proj_sol_z, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &L2_proj_mvelo);
  VecSetSizes(L2_proj_mvelo, PETSC_DECIDE, n_L2proj_totalrow);
  VecSetFromOptions(L2_proj_mvelo);
  VecSet(L2_proj_mvelo, 0.0);
  VecSetOption(L2_proj_mvelo, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  Assem_L2_proj_mat(anci.A_fixed_elementv, anci.A_elements, anci.A_quad_s, anci.A_itf_part, anci.A_SI_sol);
}

// EOF
