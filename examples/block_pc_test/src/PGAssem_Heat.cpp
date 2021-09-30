#include "PGAssem_Heat.hpp"

PGAssem_Heat::PGAssem_Heat(
    IPLocAssem * const &locassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &part_nbc,
    const ALocal_EBC * const &part_ebc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat( locassem_ptr->get_dof_mat() ),
  num_ebc( part_ebc->get_num_ebc() )
{
  SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),
      "PGAssem_HEAT::dof_sol != locassem_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != part_nbc->get_dofMat(),
      "PGAssem_HEAT::dof_mat != part_nbc->get_dofMat(). \n");

  if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  else snLocBas = 0;

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_Heat, snLocBas has to be uniform. \n");
  }

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlgn       = pnode_ptr->get_nlocghonode();
  const int nlocrow    = dof_mat * nlocalnode;

  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];
  const int empirical_neibor_number = in_nz_estimate;

  SYS_T::commPrint("     Empirical nonzero estimate: %d \n", empirical_neibor_number);

  Get_dnz_onz( nlocalnode,  empirical_neibor_number, part_nbc, dnnz, onnz );

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, dnnz, 0, onnz, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  delete [] dnnz; dnnz = nullptr;
  delete [] onnz; onnz = nullptr;

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  row_index = new PetscInt [nLocBas * dof_mat];

  array_a = new double [nlgn * dof_sol];
  array_b = new double [nlgn * dof_sol];

  local_a = new double [nLocBas * dof_sol];
  local_b = new double [nLocBas * dof_sol];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];

  if(num_ebc > 0)
  {
    LSIEN = new int [snLocBas];
    local_as = new double [dof_sol * snLocBas];
    local_bs = new double [dof_sol * snLocBas];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index = new PetscInt [dof_mat * snLocBas];
 
    VecDuplicateVecs(G, num_ebc, &intNA); 
    VecDuplicate(G, &b);
    VecDuplicate(G, &diag);
    KSPCreate(PETSC_COMM_WORLD, &ksp_K);
    KSPSetOptionsPrefix( ksp_K, "shell_" );
    KSPSetFromOptions( ksp_K ); 
  }

  for(int ii=0; ii<nlgn*dof_sol; ++ii) array_b[ii] = 0.0;

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr,
      elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ebc );

  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}


PGAssem_Heat::~PGAssem_Heat()
{
  VecDestroy(&G);
  MatDestroy(&K);
  delete [] row_index; row_index = nullptr;
  delete [] array_a;     array_a = nullptr;
  delete [] array_b;     array_b = nullptr;
  delete [] local_a;     local_a = nullptr;
  delete [] local_b;     local_b = nullptr;
  delete [] IEN_e;       IEN_e = nullptr;
  delete [] ectrl_x;     ectrl_x = nullptr;
  delete [] ectrl_y;     ectrl_y = nullptr;
  delete [] ectrl_z;     ectrl_z = nullptr;

  if(num_ebc > 0)
  {
    delete [] LSIEN; LSIEN = nullptr;
    delete [] local_as; local_as = nullptr;
    delete [] local_bs; local_bs = nullptr;
    delete [] sctrl_x; sctrl_x = nullptr;
    delete [] sctrl_y; sctrl_y = nullptr;
    delete [] sctrl_z; sctrl_z = nullptr;
    delete [] srow_index; srow_index = nullptr;
  
    VecDestroyVecs(num_ebc, &intNA);
    VecDestroy(&b);
    VecDestroy(&diag);
    KSPDestroy(&ksp_K);
  }
}


void PGAssem_Heat::Get_dnz_onz( const int &nlocnode,
    const int &empirical_neighbor_node_number,
    const ALocal_NodalBC * const &nbc_ptr,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nzbase = dof_mat * empirical_neighbor_node_number;

  Vec vdnz, vonz;
  VecCreateMPI(PETSC_COMM_WORLD, dof_mat * nlocnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof_mat * nlocnode, PETSC_DETERMINE, &vonz);
  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  VecSetOption(vdnz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(vonz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  for(int ii=0; ii<nlocnode; ++ii)
  {
    for(int mm=0; mm<dof_mat; ++mm)
    {
      const int row = nbc_ptr->get_LID(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  for(int mm=0; mm<dof_mat; ++mm)
  {
    // Check the master nodes for each d.o.f.
    const int num_master = nbc_ptr->get_Num_LPM(mm);
    for(int ii=0; ii<num_master; ++ii)
    {
      const int row = nbc_ptr->get_LocalMaster(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  VecAssemblyBegin(vdnz); VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz); VecAssemblyEnd(vonz);


  for(int mm=0; mm<dof_mat; ++mm)
  {
    // Check the Dirichlet nodes for each d.o.f.
    const int num_dir = nbc_ptr->get_Num_LD(mm);
    for(int ii=0; ii<num_dir; ++ii)
    {
      const int row = nbc_ptr->get_LDN(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
      VecSetValue(vonz, row, 0.0, INSERT_VALUES);
    }

    // Check the slave nodes for each d.o.f.
    const int num_slave = nbc_ptr->get_Num_LPS(mm);
    for(int ii=0; ii<num_slave; ++ii)
    {
      const int row = nbc_ptr->get_LPSN(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
      VecSetValue(vonz, row, 2.0, INSERT_VALUES);
    }
  }

  VecAssemblyBegin(vdnz); VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz); VecAssemblyEnd(vonz);

  PetscInt mat_length;
  VecGetSize(vdnz, &mat_length);

  const int max_dnz = dof_mat * nlocnode;
  const int max_onz = mat_length - dof_mat * nlocnode;

  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof_mat*nlocnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    if(dnz[ii] > max_dnz) dnz[ii] = max_dnz;
  }
  VecRestoreArray(vdnz, &array_d);

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof_mat*nlocnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);
    if(onz[ii] > max_onz) onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}


void PGAssem_Heat::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  lassem_ptr->Assem_Estimate();

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

  // Create a temporary zero solution vector to feed Natbc_Resis_KG
  //PDNSolution * temp = new PDNSolution_NS( node_ptr, 0, false );

  // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  //NatBC_Resis_KG(0.1, temp, temp, lassem_ptr, elements, quad_s, node_ptr,
  //    nbc_part, ebc_part, gbc );

  //delete temp;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

}


void PGAssem_Heat::Assem_tangent_residual(
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
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
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

  // Resistance type boundary condition
  //NatBC_Resis_KG( dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, node_ptr, nbc_part, ebc_part, gbc );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Heat::EssBC_KG( const ALocal_NodalBC * const &nbc_part, 
    const int &field )
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


void PGAssem_Heat::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const int &in_loc_dof,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
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

      VecSetValues(G, in_loc_dof, srow_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
}


void PGAssem_Heat::Assem_intNA(
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  PetscScalar * Res = new PetscScalar [snLocBas * 1];
  PetscInt * srow_idx = new PetscInt [snLocBas * 1];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    VecSet( intNA[ebc_id], 0.0 ); // clean the vector a^k

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
      
      lassem_ptr->Assem_Residual_EBC_Resistance( ebc_id, 1.0,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[ii] = lassem_ptr->Residual[ii];

        srow_idx[ii] = nbc_part->get_LID(0, LSIEN[ii]);
      }

      VecSetValues( intNA[ebc_id], snLocBas, srow_idx, Res, ADD_VALUES );
    }
    
    VecAssemblyBegin(intNA[ebc_id]);
    VecAssemblyEnd(intNA[ebc_id]);
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
}

// EOF
