#include "PGAssem_Wall_Prestress.hpp"

PGAssem_Wall_Prestress::PGAssem_Wall_Prestress( 
    IPLocAssem_2x2Block * const &locassem_s_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_v,
    const ALocal_IEN * const &aien_p,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const ALocal_NodalBC * const &part_nbc_v,
    const ALocal_NodalBC * const &part_nbc_p,
    const ALocal_EBC * const &part_ebc,
    const int &in_nz_estimate )
: nLocBas(4), snLocBas(3),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn_v( pnode_v -> get_nlocghonode() ),
  nlgn_p( pnode_p -> get_nlocghonode() )
{
  SYS_T::print_fatal_if( nLocBas != locassem_s_ptr->get_nLocBas_0(),
      "Error: PGAssem_FSI::nLocBas does not match that in local assembly of solid.\n");

  SYS_T::print_fatal_if( snLocBas != locassem_s_ptr->get_snLocBas_0(),
      "Error: PGAssem_FSI::nLocBas does not match that in local assembly of solid.\n");

  // Make sure the data structure is compatible
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_FSI, snLocBas has to be uniform.\n");
  }

  const int nloc_v = pnode_v -> get_nlocalnode();
  const int nloc_p = pnode_p -> get_nlocalnode();

  const int nlocrow_v = 3 * nloc_v;
  const int nlocrow_p = 1 * nloc_p;
  const int nlocrow = nlocrow_v + nlocrow_p;

  SYS_T::commPrint("     Empirical nonzero estimate: %d \n", in_nz_estimate);

  // Create vector
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  // Create matrix with routh preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 4*in_nz_estimate, PETSC_NULL, 4*in_nz_estimate, PETSC_NULL, &K);

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate( alelem_ptr, locassem_s_ptr, aien_v, aien_p, part_nbc_v, part_nbc_p );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_Wall_Prestress::~PGAssem_Wall_Prestress()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_Wall_Prestress::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const ALocal_NodalBC * const &nbc_v,
    const ALocal_NodalBC * const &nbc_p )
{
  const int nElem = alelem_ptr->get_nlocalele();

  lassem_s_ptr->Assem_Estimate();

  PetscInt * row_id_v = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      row_id_v[3*ii  ] = nbc_v -> get_LID( 0, lien_v -> get_LIEN(ee, ii) );
      row_id_v[3*ii+1] = nbc_v -> get_LID( 1, lien_v -> get_LIEN(ee, ii) );
      row_id_v[3*ii+2] = nbc_v -> get_LID( 2, lien_v -> get_LIEN(ee, ii) );
    }

    for(int ii=0; ii<nLocBas; ++ii)
      row_id_p[ii] = nbc_p -> get_LID( lien_p -> get_LIEN(ee,ii) );

    if( alelem_ptr->get_elem_tag(ee) == 1 )
    {
      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K, 3*nLocBas, row_id_v,   nLocBas, row_id_p, lassem_s_ptr->Tangent01, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent10, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_s_ptr->Tangent11, ADD_VALUES);
    }
  }

  delete [] row_id_v; row_id_v = nullptr; delete [] row_id_p; row_id_p = nullptr;

  EssBC_KG( nbc_v, nbc_p );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

void PGAssem_Wall_Prestress::EssBC_KG( const ALocal_NodalBC * const &nbc_v, 
    const ALocal_NodalBC * const &nbc_p )
{
  // For three velocity fields
  for(int field=0; field<3; ++field)
  {
    const int local_dir = nbc_v -> get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_v -> get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
        MatSetValue(K, row, row, 1.0, ADD_VALUES);
      }
    }

    const int local_sla = nbc_v -> get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_v -> get_LPSN(field, ii);
        const int col = nbc_v -> get_LPMN(field, ii);
        MatSetValue(K, row, col, 1.0, ADD_VALUES);
        MatSetValue(K, row, row, -1.0, ADD_VALUES);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }

  // For pressure field
  const int local_dir = nbc_p -> get_Num_LD(0);
  if(local_dir > 0)
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      int row = nbc_p -> get_LDN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc_p -> get_Num_LPS(0);
  if(local_sla > 0)
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      int row = nbc_p -> get_LPSN(0, ii);
      int col = nbc_p -> get_LPMN(0, ii);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_Wall_Prestress::EssBC_G( const ALocal_NodalBC * const &nbc_v, 
    const ALocal_NodalBC * const &nbc_p )
{
  // For three velocity fields
  for(int field=0; field<3; ++field)
  {
    const int local_dir = nbc_v -> get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_v -> get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc_v -> get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_v -> get_LPSN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }

  // For pressure field
  const int local_dir = nbc_p -> get_Num_LD(0);
  if(local_dir > 0)
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      int row = nbc_p -> get_LDN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_p -> get_Num_LPS(0);
  if(local_sla > 0)
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      int row = nbc_p -> get_LPSN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_Wall_Prestress::NatBC_G( const double &curr_time,
    const PDNSolution * const &pres,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part )  
{
  SYS_T::print_fatal_if( num_ebc != 1, "Error: there are more than 1 ebc surfaces.\n");
  
  const int ebc_id = 0;

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  const std::vector<double> array_p = pres -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscInt * srow_index = new PetscInt [3*snLocBas];

  for(int ee=0; ee<num_sele; ++ee)
  {
    const std::vector<int> LSIEN = ebc_part -> get_SIEN(ebc_id, ee);
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    const std::vector<double> local_p = GetLocal( array_p, LSIEN, snLocBas, 1 );

    lassem_s_ptr -> Assem_Residual_Interior_Wall_EBC( curr_time, &local_p[0], element_s, sctrl_x, sctrl_y, sctrl_z, quad_s );

    for(int ii=0; ii<snLocBas; ++ii)
    {
      srow_index[3*ii+0] = nbc_v -> get_LID(0, LSIEN[ii]);
      srow_index[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
      srow_index[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
    }

    VecSetValues(G, 3*snLocBas, srow_index, lassem_s_ptr->sur_Residual0, ADD_VALUES);
  }

  delete [] srow_index; srow_index = nullptr;
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

void PGAssem_Wall_Prestress::Assem_residual(
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_v,
    const ALocal_NodalBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const Prestress_solid * const &ps_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();

  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_dot_v = dot_velo -> GetLocalArray();
  const std::vector<double> array_dot_p = dot_pres -> GetLocalArray();

  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  PetscInt * row_id_v = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN_v = lien_v -> get_LIEN( ee );
    const std::vector<int> IEN_p = lien_p -> get_LIEN( ee );

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

    const std::vector<double> local_dot_d = GetLocal( array_dot_d, IEN_v, nLocBas, 3 );
    const std::vector<double> local_dot_v = GetLocal( array_dot_v, IEN_v, nLocBas, 3 );
    const std::vector<double> local_dot_p = GetLocal( array_dot_p, IEN_p, nLocBas, 1 );

    const std::vector<double> local_d = GetLocal( array_d, IEN_v, nLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, IEN_v, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN_p, nLocBas, 1 );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      row_id_v[3*ii  ] = nbc_v -> get_LID(0, IEN_v[ii]);
      row_id_v[3*ii+1] = nbc_v -> get_LID(1, IEN_v[ii]);
      row_id_v[3*ii+2] = nbc_v -> get_LID(2, IEN_v[ii]);
    }

    for(int ii=0; ii<nLocBas; ++ii)
      row_id_p[ii] = nbc_p -> get_LID( IEN_p[ii] );

    if( alelem_ptr->get_elem_tag(ee) == 1 )
    {
      const std::vector<double> quaprestress = ps_ptr->get_prestress( ee );

      lassem_s_ptr -> Assem_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0], 
          &local_dot_p[0], &local_d[0], &local_v[0], &local_p[0], elementv, 
          ectrl_x, ectrl_y, ectrl_z, &quaprestress[0], quad_v );

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_s_ptr->Residual0, ADD_VALUES);
      VecSetValues(G,   nLocBas, row_id_p, lassem_s_ptr->Residual1, ADD_VALUES);
    }
  }

  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  ectrl_x = nullptr; ectrl_y = nullptr; ectrl_z = nullptr;
  delete [] row_id_v; delete [] row_id_p;
  row_id_v = nullptr; row_id_p = nullptr;

  NatBC_G( curr_time, pres, lassem_s_ptr, elements, quad_s, nbc_v, ebc_part );
  
  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_G( nbc_v, nbc_p );

  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_Wall_Prestress::Assem_tangent_residual(
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_v,
    const ALocal_NodalBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const Prestress_solid * const &ps_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();

  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_dot_v = dot_velo -> GetLocalArray();
  const std::vector<double> array_dot_p = dot_pres -> GetLocalArray();

  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  PetscInt * row_id_v = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN_v = lien_v -> get_LIEN( ee );
    const std::vector<int> IEN_p = lien_p -> get_LIEN( ee );

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

    const std::vector<double> local_dot_d = GetLocal( array_dot_d, IEN_v, nLocBas, 3 );
    const std::vector<double> local_dot_v = GetLocal( array_dot_v, IEN_v, nLocBas, 3 );
    const std::vector<double> local_dot_p = GetLocal( array_dot_p, IEN_p, nLocBas, 1 );

    const std::vector<double> local_d = GetLocal( array_d, IEN_v, nLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, IEN_v, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN_p, nLocBas, 1 );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      row_id_v[3*ii  ] = nbc_v -> get_LID(0, IEN_v[ii]);
      row_id_v[3*ii+1] = nbc_v -> get_LID(1, IEN_v[ii]);
      row_id_v[3*ii+2] = nbc_v -> get_LID(2, IEN_v[ii]);
    }

    for(int ii=0; ii<nLocBas; ++ii)
      row_id_p[ii] = nbc_p -> get_LID( IEN_p[ii] );


    if( alelem_ptr->get_elem_tag(ee) == 1 )
    {
      const std::vector<double> quaprestress = ps_ptr->get_prestress( ee );

      lassem_s_ptr -> Assem_Tangent_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0],
          &local_dot_p[0], &local_d[0], &local_v[0], &local_p[0], elementv,
          ectrl_x, ectrl_y, ectrl_z, &quaprestress[0], quad_v );

      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K, 3*nLocBas, row_id_v,   nLocBas, row_id_p, lassem_s_ptr->Tangent01, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent10, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_s_ptr->Tangent11, ADD_VALUES);

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_s_ptr->Residual0, ADD_VALUES);
      VecSetValues(G,   nLocBas, row_id_p, lassem_s_ptr->Residual1, ADD_VALUES);
    }
  }

  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  ectrl_x = nullptr; ectrl_y = nullptr; ectrl_z = nullptr;
  delete [] row_id_v; delete [] row_id_p;
  row_id_v = nullptr; row_id_p = nullptr;

  NatBC_G( curr_time, pres, lassem_s_ptr, elements, quad_s, nbc_v, ebc_part );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_KG( nbc_v, nbc_p );

  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

// EOF
