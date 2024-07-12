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
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat( locassem_ptr->get_dof_mat() ),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn( pnode_ptr->get_nlocghonode() ),
  snLocBas( 0 ) 
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

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, 
      elements, elementvs, elementvs_rotated, quads, free_quad, aien_ptr, pnode_ptr, part_nbc, part_ebc, part_itf, gbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation
 
  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);

  // Allocate the vector Disp
  const int nlocrow_disp = 3 * pnode_ptr->get_nlocalnode();
  VecCreate(PETSC_COMM_WORLD, &Disp);
  VecSetSizes(Disp, nlocrow_disp, PETSC_DECIDE);

  VecSetFromOptions(Disp);
  VecSet(Disp, 0.0);
  VecSetOption(Disp, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

PGAssem_NS_FEM::~PGAssem_NS_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
  VecDestroy(&Disp);
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

  Interface_KG(0.0, 0.1, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part);

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
    const ALocal_Interface * const &itf_part )
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
  Weak_EssBC_G(0, 0, sol_a, alelem_ptr, lassem_ptr, elementvs, quad_s,
    lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // Surface integral from Nitsche method
  Interface_G(0, 0, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part);

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
    const ALocal_Interface * const &itf_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
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
  PetscInt * row_disp_index = new PetscInt [nLocBas * 3];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
    //     elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
    }

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<3; ++mm)
        row_disp_index[3*ii+mm] = 3 * nbc_part -> get_LID(0, IEN_e[ii]) + mm;
    }    
    
    if( alelem_ptr->get_elem_rotated(ee) == 0 )
    {
      lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }
    else
    {
      lassem_ptr->Assem_Residual_Rotated(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    
    VecSetValues(Disp, 3 * nLocBas, row_disp_index, lassem_ptr->disp_mesh, INSERT_VALUES);
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
  delete [] row_disp_index; row_disp_index = nullptr;
   
  // Backflow stabilization residual contribution
  BackFlow_G( sol_a, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_G( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
  //    nbc_part, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_G(curr_time, dt, sol_b, alelem_ptr, lassem_ptr, elementvs, quad_s,
      lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // For Poiseuille flow
  NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // Surface integral from Nitsche method
  Interface_G(curr_time, dt, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_G( nbc_part, ii );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  VecAssemblyBegin(Disp);
  VecAssemblyEnd(Disp);
}

void PGAssem_NS_FEM::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
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
    const ALocal_Interface * const &itf_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
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
  PetscInt * row_disp_index = new PetscInt [nLocBas * 3];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
    //     elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_e[ii])+mm;
    }

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<3; ++mm)
        row_disp_index[3*ii+mm] = 3 * nbc_part -> get_LID(0, IEN_e[ii]) + mm;
    }

    if( alelem_ptr->get_elem_rotated(ee) == 0 )
    {
      lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }
    else
    {
      lassem_ptr->Assem_Tangent_Residual_Rotated(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }
      
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
      lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);

    VecSetValues(Disp, 3 * nLocBas, row_disp_index, lassem_ptr->disp_mesh, INSERT_VALUES);
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
  delete [] row_disp_index; row_disp_index = nullptr;

  // Backflow stabilization residual & tangent contribution
  BackFlow_KG( dt, sol_a, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
  //    nbc_part, ebc_part, gbc );

  // Weakly enforced no-slip boundary condition
  // If wall_model_type = 0, it will do nothing.
  Weak_EssBC_KG(curr_time, dt, sol_b, alelem_ptr, lassem_ptr, elementvs, quad_s,
    lien_ptr, fnode_ptr, nbc_part, wbc_part);

  // For Poiseuille flow
  NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  // Surface integral from Nitsche method (Only assemble G at present)
  Interface_KG(curr_time, dt, lassem_ptr, elementvs, elementvs_rotated, elements, quad_s, free_quad, itf_part);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  VecAssemblyBegin(Disp);
  VecAssemblyEnd(Disp);
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

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }

  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

void PGAssem_NS_FEM::BackFlow_G( 
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_as = new double [dof_sol * snLocBas];
  double * local_bs = new double [dof_sol * snLocBas];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      lassem_ptr->Assem_Residual_BackFlowStab( local_as, local_bs,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
    }
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

void PGAssem_NS_FEM::BackFlow_KG( const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_as = new double [dof_sol * snLocBas];
  double * local_bs = new double [dof_sol * snLocBas];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      lassem_ptr->Assem_Tangent_Residual_BackFlowStab( dt, local_as, local_bs,
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

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
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
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );

  const int num_wele {wbc_part->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc_part->get_part_vol_ele_id(ee)};

    lien_ptr->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc_part->get_ele_face_id(ee)};

    if(alelem_ptr->get_elem_rotated(ee) == 0)
    {
      lassem_ptr->Assem_Tangent_Residual_Weak(curr_time, dt, local_b, element_vs,
        ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);
    }
    else
    {
      lassem_ptr->Assem_Tangent_Residual_Weak_Rotated(curr_time, dt, local_b, element_vs,
        ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);
    }

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_v[ii]) + mm;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index, lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

void PGAssem_NS_FEM::Weak_EssBC_G(
    const double &curr_time, const double &dt,
    const PDNSolution * const &sol,
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
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_v = new int [nLocBas];
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol->GetLocalArray( array_b );

  const int num_wele {wbc_part->get_num_ele()};

  // If wall_model_type = 0, num_wele will be 0 and this loop will be skipped.
  for(int ee{0}; ee < num_wele; ++ee)
  {
    const int local_ee_index {wbc_part->get_part_vol_ele_id(ee)};

    lien_ptr->get_LIEN(local_ee_index, IEN_v);
    GetLocal(array_b, IEN_v, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_v, ctrl_x, ctrl_y, ctrl_z);

    const int face_id {wbc_part->get_ele_face_id(ee)};
    
    if(alelem_ptr->get_elem_rotated(ee) == 0)
    {
      lassem_ptr->Assem_Residual_Weak(curr_time, dt, local_b, element_vs,
        ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);
    }
    else
    {
      lassem_ptr->Assem_Residual_Weak_Rotated(curr_time, dt, local_b, element_vs,
        ctrl_x, ctrl_y, ctrl_z, quad_s, face_id);
    }

    for(int ii{0}; ii < nLocBas; ++ii)
    {
      for(int mm{0}; mm < dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_v[ii]) + mm;
    }

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_b; array_b = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_v; IEN_v = nullptr;
  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
}

void PGAssem_NS_FEM::Interface_KG(
  const double &curr_time, const double &dt,
  IPLocAssem * const &lassem_ptr,
  FEAElement * const &fixed_elementv,
  FEAElement * const &rotated_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  IQuadPts * const &free_quad,
  const ALocal_Interface * const &itf_part )
{
  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  int * fixed_local_ien = new int [nLocBas];
  double * fixed_local_sol = new double [nLocBas * dof_sol];

  int * rotated_local_ien = new int [nLocBas];
  double * rotated_local_sol = new double [nLocBas * dof_sol];

  PetscInt * fixed_row_index = new PetscInt [nLocBas * dof_mat];
  PetscInt * rotated_row_index = new PetscInt [nLocBas * dof_mat];

  const int num_itf {itf_part->get_num_itf()};

  int ele_tag {-1};
  int rotated_ee {-1};
  std::vector<double> rotated_xi {1.0 / 3.0, 1.0 / 3.0};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    // SYS_T::commPrint("itf_id = %d\n", itf_id);
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      // SYS_T::commPrint("  fixed_ee = %d\n", ee);
      // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_nqp {quad_s->get_num_quadPts()};

      std::vector<double> R(nLocBas, 0.0);

      // Get the local ien and local sol of this fixed element
      itf_part->get_fixed_local(itf_id, ee, fixed_local_ien, fixed_local_sol);

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        fixed_elementv->get_R(qua, &R[0]);

        itf_part->get_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi);

        const int rotated_face_id {itf_part->get_rotated_face_id(itf_id, ele_tag, rotated_ee)};

        itf_part->get_rotated_ele_ctrlPts(itf_id, ele_tag, rotated_ee, curr_time, ctrl_x, ctrl_y, ctrl_z);

        free_quad->set_qp(0, rotated_xi);

        rotated_elementv->buildBasis(rotated_face_id, free_quad, ctrl_x, ctrl_y, ctrl_z);

        itf_part->get_rotated_local(itf_id, ele_tag, rotated_ee, rotated_local_ien, rotated_local_sol);

        const double qw = quad_s->get_qw(qua);

        lassem_ptr->Assem_Tangent_Residual_itf(qua, qw, dt, fixed_elementv, rotated_elementv, fixed_local_sol, rotated_local_sol, ctrl_x, ctrl_y, ctrl_z);

        for(int ii{0}; ii < nLocBas; ++ii)
        {
          for(int mm{0}; mm < dof_mat; ++mm)
          {
            fixed_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_fixed_ID(itf_id, mm, fixed_local_ien[ii]) + mm;
            rotated_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_rotated_ID(itf_id, mm, rotated_local_ien[ii]) + mm;
          }
        }

        MatSetValues(K, loc_dof, fixed_row_index, loc_dof, fixed_row_index, lassem_ptr->Tangent_ss, ADD_VALUES);
        MatSetValues(K, loc_dof, rotated_row_index, loc_dof, rotated_row_index, lassem_ptr->Tangent_rr, ADD_VALUES);

        // For static problem
        MatSetValues(K, loc_dof, fixed_row_index, loc_dof, rotated_row_index, lassem_ptr->Tangent_sr, ADD_VALUES);
        MatSetValues(K, loc_dof, rotated_row_index, loc_dof, fixed_row_index, lassem_ptr->Tangent_rs, ADD_VALUES);

        VecSetValues(G, loc_dof, fixed_row_index, lassem_ptr->Residual_s, ADD_VALUES);
        VecSetValues(G, loc_dof, rotated_row_index, lassem_ptr->Residual_r, ADD_VALUES);
        
        // VecAssemblyBegin(G);
        // VecAssemblyEnd(G);
      }
    }
  }

  delete [] fixed_local_ien; fixed_local_ien = nullptr;
  delete [] fixed_local_sol; fixed_local_sol = nullptr;
  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] rotated_local_sol; rotated_local_sol = nullptr;

  delete [] fixed_row_index; fixed_row_index = nullptr;
  delete [] rotated_row_index; rotated_row_index = nullptr;

  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
}

void PGAssem_NS_FEM::Interface_G(
  const double &curr_time, const double &dt,
  IPLocAssem * const &lassem_ptr,
  FEAElement * const &fixed_elementv,
  FEAElement * const &rotated_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  IQuadPts * const &free_quad,
  const ALocal_Interface * const &itf_part )
{
  const int loc_dof {dof_mat * nLocBas};
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  int * fixed_local_ien = new int [nLocBas];
  double * fixed_local_sol = new double [nLocBas * dof_sol];

  int * rotated_local_ien = new int [nLocBas];
  double * rotated_local_sol = new double [nLocBas * dof_sol];

  PetscInt * fixed_row_index = new PetscInt [nLocBas * dof_mat];
  PetscInt * rotated_row_index = new PetscInt [nLocBas * dof_mat];

  const int num_itf {itf_part->get_num_itf()};

  int ele_tag {-1};
  int rotated_ee {-1};
  std::vector<double> rotated_xi {1.0 / 3.0, 1.0 / 3.0};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    // SYS_T::commPrint("itf_id = %d\n", itf_id);
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      // SYS_T::commPrint("  fixed_ee = %d\n", ee);
      // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_nqp {quad_s->get_num_quadPts()};

      std::vector<double> R(nLocBas, 0.0);

      // Get the local ien and local sol of this fixed element
      itf_part->get_fixed_local(itf_id, ee, fixed_local_ien, fixed_local_sol);

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        fixed_elementv->get_R(qua, &R[0]);

        itf_part->get_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi);

        const int rotated_face_id {itf_part->get_rotated_face_id(itf_id, ele_tag, rotated_ee)};

        itf_part->get_rotated_ele_ctrlPts(itf_id, ele_tag, rotated_ee, curr_time, ctrl_x, ctrl_y, ctrl_z);

        free_quad->set_qp(0, rotated_xi);

        rotated_elementv->buildBasis(rotated_face_id, free_quad, ctrl_x, ctrl_y, ctrl_z);

        itf_part->get_rotated_local(itf_id, ele_tag, rotated_ee, rotated_local_ien, rotated_local_sol);

        const double qw = quad_s->get_qw(qua);

        lassem_ptr->Assem_Residual_itf(qua, qw, dt, fixed_elementv, rotated_elementv, fixed_local_sol, rotated_local_sol, ctrl_x, ctrl_y, ctrl_z);

        for(int ii{0}; ii < nLocBas; ++ii)
        {
          for(int mm{0}; mm < dof_mat; ++mm)
          {
            fixed_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_fixed_ID(itf_id, mm, fixed_local_ien[ii]) + mm;
            rotated_row_index[dof_mat * ii + mm] = dof_mat * itf_part->get_rotated_ID(itf_id, mm, rotated_local_ien[ii]) + mm;
          }
        }

        VecSetValues(G, loc_dof, fixed_row_index, lassem_ptr->Residual_s, ADD_VALUES);
        VecSetValues(G, loc_dof, rotated_row_index, lassem_ptr->Residual_r, ADD_VALUES);
        
        // VecAssemblyBegin(G);
        // VecAssemblyEnd(G);
      }
    }
  }

  delete [] fixed_local_ien; fixed_local_ien = nullptr;
  delete [] fixed_local_sol; fixed_local_sol = nullptr;
  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] rotated_local_sol; rotated_local_sol = nullptr;

  delete [] fixed_row_index; fixed_row_index = nullptr;
  delete [] rotated_row_index; rotated_row_index = nullptr;

  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
}

void PGAssem_NS_FEM::search_all_opposite_point(
  const double &curr_time,
  FEAElement * const &fixed_elementv,
  FEAElement * const &rotated_elementv,
  FEAElement * const &elements,
  const IQuadPts * const &quad_s,
  IQuadPts * const &free_quad,
  ALocal_Interface * const &itf_part )
{
  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];

  int * fixed_local_ien = new int [nLocBas];
  double * fixed_local_sol = new double [nLocBas * dof_sol];

  int * rotated_local_ien = new int [nLocBas];
  double * rotated_local_sol = new double [nLocBas * dof_sol];

  const int num_itf {itf_part->get_num_itf()};

  for(int itf_id{0}; itf_id<num_itf; ++itf_id)
  {
    SYS_T::commPrint("itf_id = %d\n", itf_id);
    const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

    for(int ee{0}; ee<num_fixed_elem; ++ee)
    {
      // SYS_T::commPrint("  fixed_ee = %d\n", ee);
      // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

      int ele_tag {itf_part->get_fixed_ele_tag(itf_id, ee)};

      fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

      const int fixed_face_nqp {quad_s->get_num_quadPts()};

      std::vector<double> R(nLocBas, 0.0);

      for(int qua{0}; qua<fixed_face_nqp; ++qua)
      {
        fixed_elementv->get_R(qua, &R[0]);

        // The xyz-coordinates of the quadrature point
        Vector_3 coor(0.0, 0.0, 0.0);
        for(int ii{0}; ii<nLocBas; ++ii)
        {
          coor.x() += ctrl_x[ii] * R[ii];
          coor.y() += ctrl_y[ii] * R[ii];
          coor.z() += ctrl_z[ii] * R[ii];
        }

        // SYS_T::commPrint("    point %d:\n", qua);

        int rotated_ee {0};
        search_opposite_point(curr_time, coor, itf_part, itf_id, rotated_elementv, elements, ele_tag, rotated_ee, free_quad);

        std::vector<double> rotated_xi = {free_quad->get_qp(0, 0), free_quad->get_qp(0, 1)};

        itf_part->set_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi);
      }
    }
  }

  delete [] fixed_local_ien; fixed_local_ien = nullptr;
  delete [] fixed_local_sol; fixed_local_sol = nullptr;
  delete [] rotated_local_ien; rotated_local_ien = nullptr;
  delete [] rotated_local_sol; rotated_local_sol = nullptr;

  delete [] ctrl_x; ctrl_x = nullptr;
  delete [] ctrl_y; ctrl_y = nullptr;
  delete [] ctrl_z; ctrl_z = nullptr;
}

void PGAssem_NS_FEM::search_opposite_point(
  const double &curr_time,
  const Vector_3 &fixed_pt,
  const ALocal_Interface * const &itf_part,
  const int &itf_id,
  FEAElement * rotated_elementv,
  FEAElement * elements,
  int &tag,
  int &rotated_ee,
  IQuadPts * const &rotated_xi )
  {
    bool is_found = false;

    double * volctrl_x = new double [nLocBas];
    double * volctrl_y = new double [nLocBas];
    double * volctrl_z = new double [nLocBas];

    const int snlocbas = elements->get_nLocBas();

    std::vector<double> facectrl_x(snlocbas, 0.0);
    std::vector<double> facectrl_y(snlocbas, 0.0);
    std::vector<double> facectrl_z(snlocbas, 0.0);

    int rotated_face_id = -1;
    int rotated_tag = tag;
    int num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
    // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

    for(int ee{0}; ee<num_rotated_ele; ++ee)
    {
      // SYS_T::commPrint("    search rotated ee = %d\n", ee);
      itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
      
      rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

      rotated_elementv->get_face_ctrlPts(rotated_face_id,
        volctrl_x, volctrl_y, volctrl_z,
        facectrl_x, facectrl_y, facectrl_z);

      rotated_xi->reset();
      is_found = FE_T::search_closest_point(fixed_pt, elements,
        &facectrl_x[0], &facectrl_y[0], &facectrl_z[0], rotated_xi);

      if(is_found)
      {
        rotated_ee = ee;
        // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
        break;
      }
    }

    // Second try
    if(is_found == false && tag != 0)
    {
      rotated_tag = tag - 1;

      num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee{0}; ee<num_rotated_ele; ++ee)
      {
        // SYS_T::commPrint("    search rotated ee = %d\n", ee);
        itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

        rotated_elementv->get_face_ctrlPts(rotated_face_id,
          volctrl_x, volctrl_y, volctrl_z,
          facectrl_x, facectrl_y, facectrl_z);

        rotated_xi->reset();
        is_found = FE_T::search_closest_point(fixed_pt, elements,
          &facectrl_x[0], &facectrl_y[0], &facectrl_z[0], rotated_xi);

        if(is_found)
        {
          rotated_ee = ee;
          // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
          break;
        }
      }
    }

    // Third try
    if(is_found == false && tag != itf_part->get_num_tag(itf_id) - 1)
    {
      rotated_tag = tag + 1;

      num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee{0}; ee<num_rotated_ele; ++ee)
      {
        // SYS_T::commPrint("    search rotated ee = %d\n", ee);
        itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

        rotated_elementv->get_face_ctrlPts(rotated_face_id,
          volctrl_x, volctrl_y, volctrl_z,
          facectrl_x, facectrl_y, facectrl_z);

        rotated_xi->reset();
        is_found = FE_T::search_closest_point(fixed_pt, elements,
          &facectrl_x[0], &facectrl_y[0], &facectrl_z[0], rotated_xi);

        if(is_found)
        {
          rotated_ee = ee;
          // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
          break;
        }
      }
    }

    delete [] volctrl_x; volctrl_x = nullptr;
    delete [] volctrl_y; volctrl_y = nullptr;
    delete [] volctrl_z; volctrl_z = nullptr;

    SYS_T::print_fatal_if(is_found == false,
      "Error, PGAssem_NS_GEM::search_opposite_point: cannot find opposite point.\n");

    tag = rotated_tag;
  }

// EOF
