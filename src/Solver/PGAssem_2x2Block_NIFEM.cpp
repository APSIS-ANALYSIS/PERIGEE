#include "PGAssem_2x2Block_NIFEM.hpp"

PGAssem_2x2Block_NIFEM::PGAssem_2x2Block_NIFEM(
    const int &in_nlocal_elem,
    const int &in_nlocbas_0, 
    const int &in_nlocbas_1,
    const int &in_nlocbas_geo,
    const APart_Node * const &pnode_0_ptr,
    const APart_Node * const &pnode_1_ptr,
    IPLocAssem_2x2Block * const &locassem_ptr,
    ALocal_IEN const * const &aien_0_ptr,
    ALocal_IEN const * const &aien_1_ptr,
    ALocal_IEN const * const &aien_geo_ptr,
    ALocal_NodalBC_2x2Block const * const &part_nbc,
    ALocal_EBC const * const &part_ebc )
: nlocElem( in_nlocal_elem ),
  nlocnode_0( pnode_0_ptr -> get_nlocalnode() ),
  nlgn_0( pnode_0_ptr -> get_nlocghonode() ),
  nlocnode_1( pnode_1_ptr -> get_nlocalnode() ),
  nlgn_1( pnode_1_ptr -> get_nlocghonode() ),
  dof(       locassem_ptr -> get_dof() ),
  dof_mat(   locassem_ptr -> get_dof_mat() ),
  dof_mat_0( locassem_ptr -> get_dof_mat_0() ),
  dof_mat_1( locassem_ptr -> get_dof_mat_1() ),
  nLocBas_0( in_nlocbas_0 ),
  nLocBas_1( in_nlocbas_1 ),
  nLocBas_geo( in_nlocbas_geo ),
  num_ebc( part_ebc -> get_num_ebc() ),
  loc_dof_0( dof_mat_0 * nLocBas_0 ),
  loc_dof_1( dof_mat_1 * nLocBas_1 ),
  snLocBas( part_ebc -> get_cell_nLocBas(0) )
{
  // Format check
  // lassem_ptr->dof = pnode_ptr -> dof
  if( dof != pnode_0_ptr->get_dof() || dof != pnode_1_ptr -> get_dof() )
  {
    SYS_T::commPrint("Warning: The dof from Model does not match the dof from preprocessor. \n");
    PetscPrintf(PETSC_COMM_WORLD, "         dof from Model is %d \n", dof);
    PetscPrintf(PETSC_COMM_WORLD, "         dof from preprocessor is %d and %d\n", pnode_0_ptr->get_dof(), pnode_1_ptr->get_dof() );
  }

  SYS_T::print_fatal_if( dof_mat != 4, "Error: we require dof_mat to be 4! \n" );

  SYS_T::print_fatal_if( dof_mat_0 != 1, "Error: we require dof_mat_0 to be 1!\n");

  SYS_T::print_fatal_if( dof_mat_1 != 3, "Error: we require dof_mat_1 to be 3!\n");

  // Make sure there is no essential for pressure
  SYS_T::print_fatal_if( part_nbc->get_Num_LD_0() !=0, "Error: there should be no essential BC nodes for pressure. \n");

  // Make sure part_ebc and the local assembly ebc match
  SYS_T::print_fatal_if(num_ebc != locassem_ptr->get_num_ebc_fun(), "Error: The number of ebc does not match with the number of functions implemented in the local assembly routine. \n");

  // Make sure the nLocBas and LIEN match
  SYS_T::print_fatal_if(nLocBas_0 != aien_0_ptr->get_stride(), "Error: The nLocBas_0 does not match between PGAssem_2x2Block_NIFEM and ALocal_IEN. \n");

  SYS_T::print_fatal_if(nLocBas_1 != aien_1_ptr->get_stride(), "Error: The nLocBas_0 does not match between PGAssem_2x2Block_NIFEM and ALocal_IEN. \n");

  SYS_T::print_fatal_if(nLocBas_geo != aien_geo_ptr->get_stride(), "Error: The nLocBas_0 does not match between PGAssem_2x2Block_NIFEM and ALocal_IEN. \n");

  // Assign the number of local basis functions on surface element
  if(num_ebc == 0) snLocBas = 0;

  // Make sure the surfaces are of the same type of element
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id), "Error: in PGAssem_2x2Block_NIFEM, snLocBas has to be uniform over all boundary surfaces. \n");

  // Auxiliary data structure for assembly
  row_index_0 = new PetscInt [nLocBas_0 * dof_mat_0];
  row_index_1 = new PetscInt [nLocBas_1 * dof_mat_1];

  array_sol_0 = new double [ nlgn_0 * dof_mat_0 ];
  array_sol_1 = new double [ nlgn_1 * dof_mat_1 ];

  array_dot_0 = new double [ nlgn_0 * dof_mat_0 ];
  array_dot_1 = new double [ nlgn_1 * dof_mat_1 ];

  local_sol_0 = new double [ nLocBas_0 * dof_mat_0 ];
  local_sol_1 = new double [ nLocBas_1 * dof_mat_1 ];

  local_dot_0 = new double [ nLocBas_0 * dof_mat_0 ];
  local_dot_1 = new double [ nLocBas_1 * dof_mat_1 ];

  IEN_0   = new int [nLocBas_0];
  IEN_1   = new int [nLocBas_1];
  IEN_geo = new int [nLocBas_geo];

  ectrl_x = new double [nLocBas_geo];
  ectrl_y = new double [nLocBas_geo];
  ectrl_z = new double [nLocBas_geo];

  if( num_ebc > 0 )
  {
    LSIEN = new int [snLocBas];
    local_sur_0 = new double [snLocBas * dof_mat_0];
    local_sur_1 = new double [snLocBas * dof_mat_1];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index_0 = new PetscInt [snLocBas * dof_mat_0];
    srow_index_1 = new PetscInt [snLocBas * dof_mat_1];
  }
  else
  {
    LSIEN = NULL;
    local_sur_0 = NULL; local_sur_1 = NULL;
    sctrl_x = NULL; sctrl_y = NULL; sctrl_z = NULL;
    srow_index_0 = NULL; srow_index_1 = NULL;
  }

  // Vector allocation
  VecCreate(PETSC_COMM_WORLD, &G_0);
  VecSetSizes(G_0, nlocnode_0, PETSC_DECIDE);
  VecSetFromOptions(G_0);
  VecSet(G_0, 0.0);
  VecSetOption(G_0, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  VecCreate(PETSC_COMM_WORLD, &G_1);
  VecSetSizes(G_1, 3*nlocnode_1, PETSC_DECIDE);
  VecSetFromOptions(G_1);
  VecSet(G_1, 0.0);
  VecSetOption(G_1, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  // Preallocate the matrices
  const int emp_neib_node = 100;

  MatCreateAIJ(PETSC_COMM_WORLD, nlocnode_0, nlocnode_0, PETSC_DETERMINE,
      PETSC_DETERMINE, emp_neib_node, NULL, emp_neib_node, NULL, &K_00);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocnode_0, 3*nlocnode_1, PETSC_DETERMINE,
      PETSC_DETERMINE, 3*emp_neib_node, NULL, 3*emp_neib_node, NULL, &K_01);

  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocnode_1, nlocnode_0, PETSC_DETERMINE,
      PETSC_DETERMINE, emp_neib_node, NULL, emp_neib_node, NULL, &K_10);

  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocnode_1, 3*nlocnode_1, PETSC_DETERMINE,
      PETSC_DETERMINE, 3*emp_neib_node, NULL, 3*emp_neib_node, NULL, &K_11);

  // Make one assembly
  Assem_nonzero_estimate( locassem_ptr, aien_0_ptr, aien_1_ptr, part_nbc );
 
  // Obtain the number of nonzero entries
  std::vector<int> Kdnz, Konz;

  PETSc_T::Get_dnz_onz(K_00, Kdnz, Konz);
  MatDestroy(&K_00);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocnode_0, nlocnode_0, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_00);

  PETSc_T::Get_dnz_onz(K_01, Kdnz, Konz);
  MatDestroy(&K_01);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocnode_0, 3*nlocnode_1, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_01);

  PETSc_T::Get_dnz_onz(K_10, Kdnz, Konz);
  MatDestroy(&K_10);
  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocnode_1, nlocnode_0, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_10);

  PETSc_T::Get_dnz_onz(K_11, Kdnz, Konz);
  MatDestroy(&K_11);
  MatCreateAIJ(PETSC_COMM_WORLD, 3*nlocnode_1, 3*nlocnode_1, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K_11);
}


PGAssem_2x2Block_NIFEM::~PGAssem_2x2Block_NIFEM()
{
  VecDestroy(&G_0); VecDestroy(&G_1);

  MatDestroy(&K_00); MatDestroy(&K_01);
  MatDestroy(&K_10); MatDestroy(&K_11);

  delete [] row_index_0; row_index_0 = NULL;
  delete [] row_index_1; row_index_1 = NULL;
  delete [] array_sol_0; array_sol_0 = NULL;
  delete [] array_sol_1; array_sol_1 = NULL;
  delete [] array_dot_0; array_dot_0 = NULL;
  delete [] array_dot_1; array_dot_1 = NULL;
  delete [] local_sol_0; local_sol_0 = NULL;
  delete [] local_sol_1; local_sol_1 = NULL;
  delete [] local_dot_0; local_dot_0 = NULL;
  delete [] local_dot_1; local_dot_1 = NULL;
  delete [] IEN_0; IEN_0 = NULL;
  delete [] IEN_1; IEN_1 = NULL;
  delete [] IEN_geo; IEN_geo = NULL;
  delete [] ectrl_x; ectrl_x = NULL;
  delete [] ectrl_y; ectrl_y = NULL;
  delete [] ectrl_z; ectrl_z = NULL;

  if( num_ebc > 0 )
  {
    delete [] LSIEN; LSIEN = NULL;
    delete [] local_sur_0; local_sur_0 = NULL;
    delete [] local_sur_1; local_sur_1 = NULL;
    delete [] sctrl_x; sctrl_x = NULL;
    delete [] sctrl_y; sctrl_y = NULL;
    delete [] sctrl_z; sctrl_z = NULL;
    delete [] srow_index_0; srow_index_0 = NULL;
    delete [] srow_index_1; srow_index_1 = NULL;
  }
}


void PGAssem_2x2Block_NIFEM::Get_res_norms( const PDNSolution * const &sol_0,
    const PDNSolution * const &sol_1,
    double &norm_0, double &norm_1 ) const
{
  Vec res, kxsol;
  VecDuplicate( G_1, &res ); VecDuplicate( G_1, &kxsol );

  // res = G_1 - K_11 sol_1 - K10 sol_0
  VecCopy(G_1, res);
  MatMult(K_11, sol_1->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  MatMult(K_10, sol_0->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  VecNorm(res, NORM_2, &norm_1);
  VecDestroy(&res); VecDestroy(&kxsol);

  // res = G_0 - K01 sol_1 - K00 sol_0
  VecDuplicate( G_0, &res ); VecDuplicate( G_0, &kxsol );
  VecCopy(G_0, res);
  MatMult(K_01, sol_1->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  MatMult(K_00, sol_0->solution, kxsol); VecAXPY(res, -1.0, kxsol);
  VecNorm(res, NORM_2, &norm_0);
  VecDestroy(&res); VecDestroy(&kxsol);
}


void PGAssem_2x2Block_NIFEM::EssBC_K11G1( 
    const ALocal_NodalBC_2x2Block * const &nbc_part )
{
  for(int field = 0; field < 3; ++field)
  {
    const int local_dir = nbc_part -> get_Num_LD_1( field );
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part -> get_LDN_1(field, ii) * 3 + field;
        VecSetValue(G_1, row, 0.0, INSERT_VALUES);
        MatSetValue(K_11, row, row, 1.0, ADD_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_NIFEM::EssBC_G1( 
    const ALocal_NodalBC_2x2Block * const &nbc_part )
{
  for(int field = 0; field < 3; ++field)
  {
    const int local_dir = nbc_part -> get_Num_LD_1( field );
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part ->get_LDN_1(field, ii) * 3 + field;
        VecSetValue(G_1, row, 0.0, INSERT_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_NIFEM::NatBC_G1( const double &curr_time, const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC_2x2Block * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double int_x, int_y, int_z;

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      ebc_part -> get_intPts_xyz(ebc_id, ee, int_x, int_y, int_z);

      GetLocal_0(array_sol_0, LSIEN, snLocBas, local_sur_0);
      GetLocal_1(array_sol_1, LSIEN, snLocBas, local_sur_1);

      lassem_ptr->Assem_Residual_EBC( ebc_id, curr_time, dt,
          int_x, int_y, int_z,
          local_sur_0, local_sur_1, local_sur_0, local_sur_1,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        const int loc_index = LSIEN[ii];
        for(int mm=0; mm<3; ++mm)
        {
          int lrow_index = nbc_part -> get_LID_1(mm, loc_index);
          srow_index_1[3*ii + mm] = 3 * lrow_index + mm;
        }
      }
      VecSetValues(G_1, 3*snLocBas, srow_index_1, lassem_ptr->Residual1, ADD_VALUES);
    }
  }
}


void PGAssem_2x2Block_NIFEM::GetLocal_0( const double * const &array,
    const int * const &IEN, const int &val_locbas,
    double * const &local_array) const
{
  for(int ii=0; ii<val_locbas; ++ii) local_array[ii] = array[IEN[ii]];
}


void PGAssem_2x2Block_NIFEM::GetLocal_1( const double * const &array,
    const int * const &IEN, const int &val_locbas,
    double * const &local_array) const
{
  for(int ii=0; ii<val_locbas; ++ii)
  {
    local_array[3*ii]   = array[3*IEN[ii]];
    local_array[3*ii+1] = array[3*IEN[ii]+1];
    local_array[3*ii+2] = array[3*IEN[ii]+2];
  }
}


void PGAssem_2x2Block_NIFEM::Assem_nonzero_estimate(
    IPLocAssem_2x2Block * const &lassem_ptr,
    const ALocal_IEN * const &lien_0_ptr,
    const ALocal_IEN * const &lien_1_ptr,
    const ALocal_NodalBC_2x2Block * const &nbc_part )
{
  lassem_ptr -> Assem_Estimate();

  for(int ee=0; ee<nlocElem; ++ee)
  {
    for(int ii=0; ii<nLocBas_0; ++ii)
    {
      const int loc_index = lien_0_ptr -> get_LIEN(ee, ii);

      row_index_0[ii] = nbc_part->get_LID_0(loc_index);
    }

    for(int ii=0; ii<nLocBas_1; ++ii)
    {
      const int loc_index = lien_1_ptr -> get_LIEN(ee, ii);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part->get_LID_1( mm, loc_index );

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


void PGAssem_2x2Block_NIFEM::Assem_tangent_residual(
    const PDNSolution * const &dot_sol_0,
    const PDNSolution * const &dot_sol_1,
    const PDNSolution * const &sol_0,
    const PDNSolution * const &sol_1,
    const double &curr_time,
    const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_0_ptr,
    const ALocal_IEN * const &lien_1_ptr,
    const ALocal_IEN * const &lien_geo_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC_2x2Block * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  dot_sol_0 -> GetLocalArray( array_dot_0 );
  dot_sol_1 -> GetLocalArray( array_dot_1 );

  sol_0 -> GetLocalArray( array_sol_0 );
  sol_1 -> GetLocalArray( array_sol_1 );

  for(int ee=0; ee<nlocElem; ++ee)
  {
    lien_0_ptr   -> get_LIEN_e(ee, IEN_0);
    lien_1_ptr   -> get_LIEN_e(ee, IEN_1);
    lien_geo_ptr -> get_LIEN_e(ee, IEN_geo);

    GetLocal_0(array_dot_0, IEN_0, nLocBas_0, local_dot_0);
    GetLocal_1(array_dot_1, IEN_1, nLocBas_1, local_dot_1);
    GetLocal_0(array_sol_0, IEN_0, nLocBas_0, local_sol_0);
    GetLocal_1(array_sol_1, IEN_1, nLocBas_1, local_sol_1);

    fnode_ptr->get_ctrlPts_xyz(nLocBas_geo, IEN_geo, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual( curr_time, dt, local_dot_0,
        local_dot_1, local_sol_0, local_sol_1, elementv,
        ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas_0; ++ii)
    {
      const int loc_index = lien_0_ptr -> get_LIEN(ee, ii);

      row_index_0[ii] = nbc_part->get_LID_0(loc_index);
    }

    for(int ii=0; ii<nLocBas_1; ++ii)
    {
      const int loc_index = lien_1_ptr -> get_LIEN(ee, ii);

      for(int mm=0; mm<3; ++mm)
      {
        int lrow_index = nbc_part->get_LID_1( mm, loc_index );

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

  NatBC_G1( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

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
