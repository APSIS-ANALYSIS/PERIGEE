#include "PGAssem_2x2Block_HED.hpp"

PGAssem_2x2Block_HED::PGAssem_2x2Block_HED(
    IPLocAssem_2x2Block * const &locassem_ptr,
    ALocal_IEN const * const &aien_ptr,
    APart_Node const * const &pnode_ptr,
    ALocal_NodalBC const * const &part_nbc,
    ALocal_EBC const * const &part_ebc )
: nlocElem( aien_ptr-> get_nlocalele() ),
  dof(       locassem_ptr -> get_dof() ),
  dof_mat(   locassem_ptr -> get_dof_mat() ),
  dof_mat_0( locassem_ptr -> get_dof_mat_0() ),
  dof_mat_1( locassem_ptr -> get_dof_mat_1() ),
  nLocBas( aien_ptr-> get_stride() ),
  num_ebc( part_ebc -> get_num_ebc() ),
  loc_dof_0( dof_mat_0 * nLocBas ),
  loc_dof_1( dof_mat_1 * nLocBas )
{
  // We pass nLocBas and nlocElem from ALocal_IEN. This is because
  // we assume equal order interpolation now. For non-isoparametric
  // element, one should be careful.
  
  // Format check
  if( dof != pnode_ptr->get_dof() )
  {
    SYS_T::commPrint("Warning: The dof from Model does not match the dof from preprocessor. \n");
    PetscPrintf(PETSC_COMM_WORLD, "         dof from Model is %d \n", dof);
    PetscPrintf(PETSC_COMM_WORLD, "         dof from preprocessor is %d \n", pnode_ptr->get_dof());
  }

  SYS_T::print_fatal_if( dof_mat != 4 , "Error: we require dof_mat to be 4! \n" );

  SYS_T::print_fatal_if( dof_mat_0 != 1, "Error: we require dof_mat_0 to be 1!\n");

  SYS_T::print_fatal_if( dof_mat_1 != 3, "Error: we require dof_mat_1 to be 3!\n");

  for(int field=0; field<4; ++field) SYS_T::print_fatal_if( part_nbc-> get_Num_LPS(field) != 0, "Error: periodic boundary condition is NOT implemented in PGAssem_2x2Block_FEM. \n");

  SYS_T::print_fatal_if( part_nbc->get_Num_LD(0) !=0, "Error: there should be no essential BC nodes for pressure. \n");

  SYS_T::print_fatal_if(num_ebc != locassem_ptr->get_num_ebc_fun(), "Error: The number of ebc does not match with the number of functions implemented in the local assembly routine. \n");

  // Assign surface element basis function number
  if(num_ebc > 0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  else snLocBas = 0;

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id), "Error: in PGAssem_2x2Block_FEM, snLocBas has to be uniform. \n");

  // Auxiliary variable for assembly
  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlgn       = pnode_ptr->get_nlocghonode();
  
  row_index_0 = new PetscInt [nLocBas * dof_mat_0];
  row_index_1 = new PetscInt [nLocBas * dof_mat_1];
  
  array_sol_d = new double [ nlgn * dof_mat_1 ];
  array_sol_p = new double [ nlgn * dof_mat_0 ];
  array_sol_v = new double [ nlgn * dof_mat_1 ];

  array_dot_d = new double [ nlgn * dof_mat_1 ];
  array_dot_p = new double [ nlgn * dof_mat_0 ];
  array_dot_v = new double [ nlgn * dof_mat_1 ];

  local_sol_d = new double [ nLocBas * dof_mat_1 ];
  local_sol_p = new double [ nLocBas * dof_mat_0 ];
  local_sol_v = new double [ nLocBas * dof_mat_1 ];

  local_dot_d = new double [ nLocBas * dof_mat_1 ];
  local_dot_p = new double [ nLocBas * dof_mat_0 ];
  local_dot_v = new double [ nLocBas * dof_mat_1 ];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];

  if( num_ebc > 0 )
  {
    LSIEN = new int [snLocBas];
    local_sur_d = new double [snLocBas * dof_mat_1];
    local_sur_p = new double [snLocBas * dof_mat_0];
    local_sur_v = new double [snLocBas * dof_mat_1];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index_0 = new PetscInt [snLocBas * dof_mat_0];
    srow_index_1 = new PetscInt [snLocBas * dof_mat_1];
  }
  else
  {
    LSIEN = NULL;
    local_sur_d = NULL; local_sur_p = NULL; local_sur_v = NULL;
    sctrl_x = NULL; sctrl_y = NULL; sctrl_z = NULL;
    srow_index_0 = NULL; srow_index_1 = NULL;
  }

  // Allocate vectors
  VecCreate(PETSC_COMM_WORLD, &G_0);
  VecSetSizes(G_0, nlocalnode, PETSC_DECIDE);
  VecSetFromOptions(G_0);
  VecSet(G_0, 0.0);
  VecSetOption(G_0, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &G_1);
  VecSetSizes(G_1, 3*nlocalnode, PETSC_DECIDE);
  VecSetFromOptions(G_1);
  VecSet(G_1, 0.0);
  VecSetOption(G_1, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  // Pre-allocate matrices
  const int emp_neib_node = 30;

  MatCreateAIJ(PETSC_COMM_WORLD, nlocalnode, nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, emp_neib_node, NULL, emp_neib_node, NULL, &K_00);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocalnode, 3*nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 3*emp_neib_node, NULL, 3*emp_neib_node, NULL, &K_01);

  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocalnode, nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, emp_neib_node, NULL, emp_neib_node, NULL, &K_10);

  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocalnode, 3*nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 3*emp_neib_node, NULL, 3*emp_neib_node, NULL, &K_11);

  // Generation of nonzero beyond the estimate will NOT throw an error.
  // The non-zero pattern will be fixed after we get correct allocation.
  Release_nonzero_err_str();

  // Non-zero estimate
  Assem_nonzero_estimate( locassem_ptr, aien_ptr, part_nbc );
  
  // Obtain nonzero pattern and re-generate the matrices
  std::vector<int> Kdnz, Konz;

  PETSc_T::Get_dnz_onz(K_00, Kdnz, Konz);
  MatDestroy(&K_00);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocalnode, nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_00);

  PETSc_T::Get_dnz_onz(K_01, Kdnz, Konz);
  MatDestroy(&K_01);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocalnode, 3*nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_01);

  PETSc_T::Get_dnz_onz(K_10, Kdnz, Konz);
  MatDestroy(&K_10);
  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocalnode, nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_10);

  PETSc_T::Get_dnz_onz(K_11, Kdnz, Konz);
  MatDestroy(&K_11);
  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocalnode, 3*nlocalnode, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_11);
}


PGAssem_2x2Block_HED::~PGAssem_2x2Block_HED()
{
  VecDestroy(&G_0); VecDestroy(&G_1);

  MatDestroy(&K_00); MatDestroy(&K_01);
  MatDestroy(&K_10); MatDestroy(&K_11);

  delete [] row_index_0; row_index_0 = NULL;
  delete [] row_index_1; row_index_1 = NULL;
  delete [] array_sol_d; array_sol_d = NULL;
  delete [] array_sol_p; array_sol_p = NULL;
  delete [] array_sol_v; array_sol_v = NULL;
  delete [] array_dot_d; array_dot_d = NULL;
  delete [] array_dot_p; array_dot_p = NULL;
  delete [] array_dot_v; array_dot_v = NULL;
  delete [] local_sol_d; local_sol_d = NULL;
  delete [] local_sol_p; local_sol_p = NULL;
  delete [] local_sol_v; local_sol_v = NULL;
  delete [] local_dot_d; local_dot_d = NULL;
  delete [] local_dot_p; local_dot_p = NULL;
  delete [] local_dot_v; local_dot_v = NULL;
  delete [] IEN_e; IEN_e = NULL;
  delete [] ectrl_x; ectrl_x = NULL;
  delete [] ectrl_y; ectrl_y = NULL;
  delete [] ectrl_z; ectrl_z = NULL;

  if( num_ebc > 0 )
  {
    delete [] LSIEN; LSIEN = NULL;
    delete [] local_sur_d; local_sur_d = NULL;
    delete [] local_sur_p; local_sur_p = NULL;
    delete [] local_sur_v; local_sur_v = NULL;
    delete [] sctrl_x; sctrl_x = NULL;
    delete [] sctrl_y; sctrl_y = NULL;
    delete [] sctrl_z; sctrl_z = NULL;
    delete [] srow_index_0; srow_index_0 = NULL;
    delete [] srow_index_1; srow_index_1 = NULL;
  }
}


void PGAssem_2x2Block_HED::Get_res_norms( 
    const PDNSolution * const &sol_0,
    const PDNSolution * const &sol_1,
    double &norm_0, double &norm_1 ) const
{
  Vec res, kxsol;
  VecDuplicate( G_1, &res ); VecDuplicate( G_1, &kxsol );

  VecCopy(G_1, res);
  MatMult(K_11, sol_1->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  MatMult(K_10, sol_0->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  VecNorm(res, NORM_2, &norm_1);
  VecDestroy(&res); VecDestroy(&kxsol);

  VecDuplicate( G_0, &res ); VecDuplicate( G_0, &kxsol );
  VecCopy(G_0, res);
  MatMult(K_01, sol_1->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  MatMult(K_00, sol_0->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  VecNorm(res, NORM_2, &norm_0);
  VecDestroy(&res); VecDestroy(&kxsol);
}


void PGAssem_2x2Block_HED::EssBC_K11G1( const ALocal_NodalBC * const &nbc_part )
{
  for(int field = 1; field <= 3; ++field)
  {
    const int local_dir = nbc_part -> get_Num_LD( field );
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = nbc_part ->get_LDN(field, ii) * 3 + field - 1;
        VecSetValue(G_1,  row, 0.0, INSERT_VALUES);
        MatSetValue(K_11, row, row, 1.0, ADD_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_HED::EssBC_G1( const ALocal_NodalBC * const &nbc_part )
{
  for(int field = 1; field <= 3; ++field)
  {
    const int local_dir = nbc_part -> get_Num_LD( field );
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = nbc_part ->get_LDN(field, ii) * 3 + field - 1;
        VecSetValue(G_1, row, 0.0, INSERT_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_HED::NatBC_G1( const double &curr_time, const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const lien_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double int_x, int_y, int_z;

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      ebc_part -> get_intPts_xyz(ebc_id, ee, int_x, int_y, int_z);

      GetLocal_0(array_sol_p, LSIEN, snLocBas, local_sur_p);
      GetLocal_1(array_sol_d, LSIEN, snLocBas, local_sur_d);
      GetLocal_1(array_sol_v, LSIEN, snLocBas, local_sur_v);

      lassem_ptr->Assem_Residual_EBC( ebc_id, curr_time, dt,
          int_x, int_y, int_z, local_sur_d, local_sur_p, local_sur_v,
          local_sur_d, local_sur_p, local_sur_v,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        int loc_index = LSIEN[ii];
        for(int mm=0; mm<3; ++mm)
        {
          int lrow_index = nbc_part -> get_LID(mm+1, loc_index);
          srow_index_1[3*ii + mm] = 3 * lrow_index + mm;
        }
      }
      VecSetValues(G_1, 3*snLocBas, srow_index_1, lassem_ptr->Residual1, ADD_VALUES);
    }
  }
}


void PGAssem_2x2Block_HED::Assem_nonzero_estimate(
    IPLocAssem_2x2Block * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_NodalBC * const &nbc_part )
{
  lassem_ptr -> Assem_Estimate();

  for(int ee=0; ee<nlocElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = lien_ptr -> get_LIEN(ee, ii);

      row_index_0[ii] = nbc_part->get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part -> get_LID( mm+1, loc_index );

        row_index_1[3*ii+mm] = 3 * lrow_index + mm;
      }
    }

    MatSetValues(K_00, loc_dof_0, row_index_0, loc_dof_0, row_index_0,
        lassem_ptr->Tangent00, ADD_VALUES);

    MatSetValues(K_01, loc_dof_0, row_index_0, loc_dof_1, row_index_1,
        lassem_ptr->Tangent01, ADD_VALUES);

    MatSetValues(K_10, loc_dof_1, row_index_1, loc_dof_0, row_index_0,
        lassem_ptr->Tangent10, ADD_VALUES);

    MatSetValues(K_11, loc_dof_1, row_index_1, loc_dof_1, row_index_1,
        lassem_ptr->Tangent11, ADD_VALUES);
  }

  EssBC_K11G1( nbc_part );

  MatAssemblyBegin(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_11, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_11, MAT_FINAL_ASSEMBLY);

  VecAssemblyBegin(G_0); VecAssemblyEnd(G_0);
  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);
}


void PGAssem_2x2Block_HED::Assem_mass_residual(
    const PDNSolution * const &sol_d,
    const PDNSolution * const &sol_p,
    const PDNSolution * const &sol_v,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  sol_d->GetLocalArray( array_sol_d );
  sol_p->GetLocalArray( array_sol_p );
  sol_v->GetLocalArray( array_sol_v );

  for(int ee=0; ee<nlocElem; ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    
    GetLocal_0(array_sol_p, IEN_e, nLocBas, local_sol_p);
    
    GetLocal_1(array_sol_d, IEN_e, nLocBas, local_sol_d);
    
    GetLocal_1(array_sol_v, IEN_e, nLocBas, local_sol_v);
    
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Mass_Residual( local_sol_d, local_sol_p,
        local_sol_v, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      int loc_index = IEN_e[ii];

      row_index_0[ii] = nbc_part->get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part->get_LID(mm+1, loc_index);
        row_index_1[3*ii+mm] = 3 * lrow_index + mm;
      }
    }

    MatSetValues(K_00, loc_dof_0, row_index_0, loc_dof_0, row_index_0,
        lassem_ptr->Tangent00, ADD_VALUES);

    MatSetValues(K_01, loc_dof_0, row_index_0, loc_dof_1, row_index_1,
        lassem_ptr->Tangent01, ADD_VALUES);

    MatSetValues(K_10, loc_dof_1, row_index_1, loc_dof_0, row_index_0,
        lassem_ptr->Tangent10, ADD_VALUES);

    MatSetValues(K_11, loc_dof_1, row_index_1, loc_dof_1, row_index_1,
        lassem_ptr->Tangent11, ADD_VALUES);

    VecSetValues(G_0, loc_dof_0, row_index_0, lassem_ptr->Residual0, ADD_VALUES);
    VecSetValues(G_1, loc_dof_1, row_index_1, lassem_ptr->Residual1, ADD_VALUES);
  }

  NatBC_G1( 0.0, 0.0, lassem_ptr, elements, quad_s, lien_ptr,
      nbc_part, ebc_part );

  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);

  EssBC_K11G1( nbc_part );

  MatAssemblyBegin(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_11, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_11, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G_0); VecAssemblyEnd(G_0);
  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);
}


void PGAssem_2x2Block_HED::Assem_residual(
    const PDNSolution * const &dot_d,
    const PDNSolution * const &dot_p,
    const PDNSolution * const &dot_v,
    const PDNSolution * const &sol_d,
    const PDNSolution * const &sol_p,
    const PDNSolution * const &sol_v,
    const double &curr_time,
    const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  dot_d -> GetLocalArray( array_dot_d );
  dot_p -> GetLocalArray( array_dot_p );
  dot_v -> GetLocalArray( array_dot_v );

  sol_d -> GetLocalArray( array_sol_d );
  sol_p -> GetLocalArray( array_sol_p );
  sol_v -> GetLocalArray( array_sol_v );

  for(int ee=0; ee<nlocElem; ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    
    GetLocal_1(array_dot_d, IEN_e, nLocBas, local_dot_d);
    GetLocal_0(array_dot_p, IEN_e, nLocBas, local_dot_p);
    GetLocal_1(array_dot_v, IEN_e, nLocBas, local_dot_v);
    
    GetLocal_1(array_sol_d, IEN_e, nLocBas, local_sol_d);
    GetLocal_0(array_sol_p, IEN_e, nLocBas, local_sol_p);
    GetLocal_1(array_sol_v, IEN_e, nLocBas, local_sol_v);
    
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Residual( curr_time, dt, local_dot_d,
        local_dot_p, local_dot_v, local_sol_d, local_sol_p,
        local_sol_v, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v ); 

    // Generate the slots for assembly. 
    for(int ii=0; ii<nLocBas; ++ii)
    {
      int loc_index = IEN_e[ii];
      row_index_0[ii] = nbc_part->get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part->get_LID(mm+1, loc_index);
        row_index_1[3*ii+mm] = 3 * lrow_index + mm;
      }
    }

    VecSetValues(G_0, loc_dof_0, row_index_0, lassem_ptr->Residual0, ADD_VALUES);
    VecSetValues(G_1, loc_dof_1, row_index_1, lassem_ptr->Residual1, ADD_VALUES);
  } 

  NatBC_G1( curr_time, dt, lassem_ptr, elements, quad_s, lien_ptr,
      nbc_part, ebc_part );

  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);

  EssBC_G1( nbc_part );

  VecAssemblyBegin(G_0); VecAssemblyEnd(G_0);
  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);
}


void PGAssem_2x2Block_HED::Assem_tangent_residual(
    const PDNSolution * const &dot_d,
    const PDNSolution * const &dot_p,
    const PDNSolution * const &dot_v,
    const PDNSolution * const &sol_d,
    const PDNSolution * const &sol_p,
    const PDNSolution * const &sol_v,
    const double &curr_time,
    const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  dot_d -> GetLocalArray( array_dot_d );
  dot_p -> GetLocalArray( array_dot_p );
  dot_v -> GetLocalArray( array_dot_v );

  sol_d -> GetLocalArray( array_sol_d );
  sol_p -> GetLocalArray( array_sol_p );
  sol_v -> GetLocalArray( array_sol_v );

  for(int ee=0; ee<nlocElem; ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    
    GetLocal_1(array_dot_d, IEN_e, nLocBas, local_dot_d);
    GetLocal_0(array_dot_p, IEN_e, nLocBas, local_dot_p);
    GetLocal_1(array_dot_v, IEN_e, nLocBas, local_dot_v);
    
    GetLocal_1(array_sol_d, IEN_e, nLocBas, local_sol_d);
    GetLocal_0(array_sol_p, IEN_e, nLocBas, local_sol_p);
    GetLocal_1(array_sol_v, IEN_e, nLocBas, local_sol_v);
    
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual( curr_time, dt, local_dot_d,
        local_dot_p, local_dot_v, local_sol_d, local_sol_p,
        local_sol_v, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v );

    // Generate the slots for assembly
    for(int ii=0; ii<nLocBas; ++ii)
    {
      int loc_index = IEN_e[ii];
      row_index_0[ii] = nbc_part->get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part->get_LID(mm+1, loc_index);
        row_index_1[3*ii+mm] = 3 * lrow_index + mm;
      }
    }

    MatSetValues(K_00, loc_dof_0, row_index_0, loc_dof_0, row_index_0,
        lassem_ptr->Tangent00, ADD_VALUES);

    MatSetValues(K_01, loc_dof_0, row_index_0, loc_dof_1, row_index_1,
        lassem_ptr->Tangent01, ADD_VALUES);

    MatSetValues(K_10, loc_dof_1, row_index_1, loc_dof_0, row_index_0,
        lassem_ptr->Tangent10, ADD_VALUES);

    MatSetValues(K_11, loc_dof_1, row_index_1, loc_dof_1, row_index_1,
        lassem_ptr->Tangent11, ADD_VALUES);

    VecSetValues(G_0, loc_dof_0, row_index_0, lassem_ptr->Residual0, ADD_VALUES);
    VecSetValues(G_1, loc_dof_1, row_index_1, lassem_ptr->Residual1, ADD_VALUES);
  }

  NatBC_G1( curr_time, dt, lassem_ptr, elements, quad_s, lien_ptr,
      nbc_part, ebc_part );

  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);

  EssBC_K11G1( nbc_part );

  MatAssemblyBegin(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_00, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_01, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_10, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K_11, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_11, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G_0); VecAssemblyEnd(G_0);
  VecAssemblyBegin(G_1); VecAssemblyEnd(G_1);
}

// EOF
