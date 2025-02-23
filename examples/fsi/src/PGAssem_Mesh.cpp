#include "PGAssem_Mesh.hpp"

PGAssem_Mesh::PGAssem_Mesh( IPLocAssem * const &locassem_ptr,
    ALocal_Elem const * const &alelem_ptr,
    ALocal_IEN const * const &aien_ptr,
    APart_Node const * const &pnode_ptr,
    ALocal_NBC const * const &part_nbc,
    ALocal_EBC const * const &part_ebc,
    const int &in_nz_estimate )
: nLocBas(  locassem_ptr->get_nLocBas() ), 
  snLocBas( locassem_ptr->get_snLocBas() ), 
  dof(3),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn( pnode_ptr->get_nlocghonode() )
{
  SYS_T::print_fatal_if(num_ebc != locassem_ptr->get_num_ebc_fun(), "Error: The number of ebc does not match with the number of functions implemented in the local assembly routine. \n");
  
  SYS_T::print_fatal_if(dof != locassem_ptr->get_dof(), "PGAssem_Mesh::dof != locassem_ptr->get_dof(). \n");

  const int nlocrow = dof * pnode_ptr->get_nlocalnode();

  // Allocate the sparse matrix K
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof*in_nz_estimate, NULL, dof*in_nz_estimate, NULL, &K);

  // Allocate the vector G
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, aien_ptr, part_nbc );
  
  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_Mesh::~PGAssem_Mesh()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_Mesh::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_NBC * const &nbc_part )
{
  const int nElem = alelem_ptr -> get_nlocalele();

  lassem_ptr->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = lien_ptr -> get_LIEN(ee, ii);

      for(int mm=0; mm<dof; ++mm)
        row_index[dof*ii+mm] = nbc_part -> get_LID(mm, loc_index);
    }

    MatSetValues(K, dof*nLocBas, row_index, dof*nLocBas, row_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  EssBC_KG( nbc_part );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Mesh::Assem_mass_residual(
    const PDNSolution * const &sol,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr -> get_nlocalele();
  const int loc_dof = dof * nLocBas;
  
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof];
  
  const std::vector<double> array = sol -> GetLocalArray();
  
  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN_e = lien_ptr -> get_LIEN(ee);
    const std::vector<double> local = GetLocal(array, IEN_e, nLocBas);
    
    fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);
    
    lassem_ptr -> Assem_Mass_Residual( &local[0], ectrl_x, ectrl_y, ectrl_z );
    
    for(int ii=0; ii<nLocBas; ++ii)
      for(int mm=0; mm<dof; ++mm)
        row_index[dof*ii+mm] = nbc_part -> get_LID(mm, IEN_e[ii]);

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index, lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;
  
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  EssBC_KG( nbc_part );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Mesh::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof];

  const std::vector<double> array_a = sol_a->GetLocalArray();
  const std::vector<double> array_b = sol_b->GetLocalArray();

  for( int ee=0; ee<nElem; ++ee )
  {
    const std::vector<int> IEN_e = lien_ptr -> get_LIEN(ee);
    const std::vector<double> local_a = GetLocal(array_a, IEN_e, nLocBas);
    const std::vector<double> local_b = GetLocal(array_b, IEN_e, nLocBas);
  
    fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);
   
    lassem_ptr->Assem_Residual(curr_time, dt, &local_a[0], &local_b[0],
		    ectrl_x, ectrl_y, ectrl_z );

    for(int ii=0; ii<nLocBas; ++ii)
      for(int mm=0; mm<dof; ++mm)
        row_index[dof*ii+mm] = nbc_part -> get_LID(mm, IEN_e[ii]);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  NatBC_G(curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  EssBC_G( nbc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Mesh::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof];

  const std::vector<double> array_a = sol_a->GetLocalArray();
  const std::vector<double> array_b = sol_b->GetLocalArray();

  for( int ee=0; ee<nElem; ++ee )
  {
   const std::vector<int> IEN_e = lien_ptr -> get_LIEN(ee);
   const std::vector<double> local_a = GetLocal(array_a, IEN_e, nLocBas);
   const std::vector<double> local_b = GetLocal(array_b, IEN_e, nLocBas);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual(curr_time, dt, &local_a[0], &local_b[0],
		    ectrl_x, ectrl_y, ectrl_z );

    for(int ii=0; ii<nLocBas; ++ii)
      for(int mm=0; mm<dof; ++mm)
        row_index[dof*ii+mm] = nbc_part -> get_LID(mm, IEN_e[ii]);

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  NatBC_G(curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_KG( nbc_part );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_Mesh::EssBC_KG( const ALocal_NBC * const &nbc_part )
{
  for(int field=0; field<dof; ++field)
  {
    // Number of local dirichlet nodes for this field
    const int local_dir = nbc_part -> get_Num_LD(field);

    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part->get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
        MatSetValue(K, row, row, 1.0, ADD_VALUES);
      }
    }

    // Number of slave nodes for this field
    const int local_sla = nbc_part->get_Num_LPS(field);
    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_part->get_LPSN(field, ii);
        const int col = nbc_part->get_LPMN(field, ii);
        MatSetValue(K, row, col, 1.0, ADD_VALUES);
        MatSetValue(K, row, row, -1.0, ADD_VALUES);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }
}

void PGAssem_Mesh::EssBC_G( const ALocal_NBC * const &nbc_part )
{
  for(int field=0; field<dof; ++field)
  {
    // Number of local dirichlet nodes for this field
    const int local_dir = nbc_part -> get_Num_LD(field);

    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part->get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }

    // Number of slave nodes for this field
    const int local_sla = nbc_part->get_Num_LPS(field);
    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_part->get_LPSN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }
}

void PGAssem_Mesh::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof * snLocBas];

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      const std::vector<int> LSIEN = ebc_part -> get_SIEN(ebc_id, ee); 
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          sctrl_x, sctrl_y, sctrl_z);

      for(int ii=0; ii<snLocBas; ++ii)
        for(int mm=0; mm<dof; ++mm)
          srow_index[dof * ii + mm] =  nbc_part -> get_LID(mm, LSIEN[ii]);

      VecSetValues(G, dof*snLocBas, srow_index, lassem_ptr->Residual, ADD_VALUES);  
    }
  }

  delete [] sctrl_x;    sctrl_x = nullptr;
  delete [] sctrl_y;    sctrl_y = nullptr;
  delete [] sctrl_z;    sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

// EOF
