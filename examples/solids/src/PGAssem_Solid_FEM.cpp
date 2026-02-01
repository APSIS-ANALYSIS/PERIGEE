#include "PGAssem_Solid_FEM.hpp"

PGAssem_Solid_FEM::PGAssem_Solid_FEM(
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    std::unique_ptr<ALocal_NBC> in_nbc,
    std::unique_ptr<ALocal_EBC> in_ebc,
    std::unique_ptr<IPLocAssem_2x2Block> in_locassem,
    const int &in_nz_estimate )
: locien( std::move(in_locien) ),
  locelem( std::move(in_locelem) ),
  fnode( std::move(in_fnode) ),
  pnode( std::move(in_pnode) ),
  nbc( std::move(in_nbc) ),
  ebc( std::move(in_ebc) ),
  locassem( std::move(in_locassem) ),
  num_ebc( ebc->get_num_ebc() ),
  nLocBas( locassem->get_nLocBas_0() ),
  snLocBas( locassem->get_snLocBas_0() ),
  dof_mat( locassem->get_dof_0() + locassem->get_dof_1() ),
  nlgn( pnode->get_nlocghonode() ),
  nqpv( locassem->get_nqpv() )
{
  SYS_T::print_fatal_if( dof_mat != nbc->get_dof_LID(),
      "Error: PGAssem_Solid_FEM, dof_mat and nbc dof mismatch.\n" );

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    SYS_T::print_fatal_if( snLocBas != ebc->get_cell_nLocBas(ebc_id),
        "Error: PGAssem_Solid_FEM, snLocBas has to be uniform.\n" );
  }

  zero_prestress.assign( nqpv * 6, 0.0 );

  const int nlocrow = dof_mat * pnode -> get_nlocalnode();

  // Create vector
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  // Create matrix with rough preallocation
#if PETSC_VERSION_LT(3,19,0)
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat*in_nz_estimate, PETSC_NULL,
      dof_mat*in_nz_estimate, PETSC_NULL, &K);
#else
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat*in_nz_estimate, PETSC_NULLPTR,
      dof_mat*in_nz_estimate, PETSC_NULLPTR, &K);
#endif

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

PGAssem_Solid_FEM::~PGAssem_Solid_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_Solid_FEM::EssBC_KG( const int &field )
{
  const int local_dir = nbc -> get_Num_LD(field);

  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc -> get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc -> get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc -> get_LPSN(field, ii) * dof_mat + field;
      const int col = nbc -> get_LPMN(field, ii) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_Solid_FEM::EssBC_G( const int &field )
{
  const int local_dir = nbc -> get_Num_LD(field);

  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc -> get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc -> get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc -> get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_Solid_FEM::Assem_nonzero_estimate()
{
  const int nElem = locelem->get_nlocalele();

  locassem->Assem_Estimate();

  PetscInt * row_id_u = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = locien -> get_LIEN(ee, ii);

      row_id_p[ii] = dof_mat * nbc -> get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
      {
        row_id_u[3*ii + mm] = dof_mat * nbc -> get_LID(mm+1, loc_index) + mm + 1;
      }
    }

    MatSetValues(K, 3*nLocBas, row_id_u, 3*nLocBas, row_id_u,
        locassem->Tangent00, ADD_VALUES);

    MatSetValues(K, 3*nLocBas, row_id_u,   nLocBas, row_id_p,
        locassem->Tangent01, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_u,
        locassem->Tangent10, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p,
        locassem->Tangent11, ADD_VALUES);
  }

  delete [] row_id_u; row_id_u = nullptr;
  delete [] row_id_p; row_id_p = nullptr;

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG(ii);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_Solid_FEM::NatBC_G( const double &curr_time, const double &dt )
{
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [3 * snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc -> get_SIEN(ebc_id, ee, LSIEN);
      ebc -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      locassem->Assem_Residual_EBC( ebc_id, curr_time, dt,
          sctrl_x, sctrl_y, sctrl_z );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<3; ++mm)
        {
          srow_index[3*ii + mm] = dof_mat * nbc -> get_LID(mm+1, LSIEN[ii]) + mm + 1;
        }
      }

      VecSetValues(G, 3*snLocBas, srow_index, locassem->sur_Residual0, ADD_VALUES);
    }
  }

  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

void PGAssem_Solid_FEM::Assem_mass_residual(
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres )
{
  const int nElem = locelem->get_nlocalele();

  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  PetscInt * row_id_u = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN = locien -> get_LIEN( ee );
    fnode -> get_ctrlPts_xyz(nLocBas, &IEN[0], ectrl_x, ectrl_y, ectrl_z);

    const std::vector<double> local_d = GetLocal( array_d, IEN, nLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, IEN, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN, nLocBas, 1 );

    locassem -> Assem_Mass_Residual( &local_d[0], &local_v[0], &local_p[0],
        ectrl_x, ectrl_y, ectrl_z, &zero_prestress[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = IEN[ii];
      row_id_p[ii] = dof_mat * nbc -> get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
        row_id_u[3*ii + mm] = dof_mat * nbc -> get_LID(mm+1, loc_index) + mm + 1;
    }

    MatSetValues(K, 3*nLocBas, row_id_u, 3*nLocBas, row_id_u,
        locassem->Tangent00, ADD_VALUES);

    MatSetValues(K, 3*nLocBas, row_id_u,   nLocBas, row_id_p,
        locassem->Tangent01, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_u,
        locassem->Tangent10, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p,
        locassem->Tangent11, ADD_VALUES);

    VecSetValues(G, 3*nLocBas, row_id_u, locassem->Residual0, ADD_VALUES);
    VecSetValues(G,   nLocBas, row_id_p, locassem->Residual1, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_id_u; row_id_u = nullptr;
  delete [] row_id_p; row_id_p = nullptr;

  NatBC_G( 0.0, 0.0 );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG(ii);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_Solid_FEM::Assem_Residual(
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres )
{
  const int nElem = locelem->get_nlocalele();

  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_dot_v = dot_velo -> GetLocalArray();
  const std::vector<double> array_dot_p = dot_pres -> GetLocalArray();

  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  PetscInt * row_id_u = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN = locien -> get_LIEN( ee );
    fnode -> get_ctrlPts_xyz(nLocBas, &IEN[0], ectrl_x, ectrl_y, ectrl_z);

    const std::vector<double> local_dot_d = GetLocal( array_dot_d, IEN, nLocBas, 3 );
    const std::vector<double> local_dot_v = GetLocal( array_dot_v, IEN, nLocBas, 3 );
    const std::vector<double> local_dot_p = GetLocal( array_dot_p, IEN, nLocBas, 1 );

    const std::vector<double> local_d = GetLocal( array_d, IEN, nLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, IEN, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN, nLocBas, 1 );

    locassem -> Assem_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0],
        &local_dot_p[0], &local_d[0], &local_v[0], &local_p[0],
        ectrl_x, ectrl_y, ectrl_z, &zero_prestress[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = IEN[ii];
      row_id_p[ii] = dof_mat * nbc -> get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
        row_id_u[3*ii + mm] = dof_mat * nbc -> get_LID(mm+1, loc_index) + mm + 1;
    }

    VecSetValues(G, 3*nLocBas, row_id_u, locassem->Residual0, ADD_VALUES);
    VecSetValues(G,   nLocBas, row_id_p, locassem->Residual1, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_id_u; row_id_u = nullptr;
  delete [] row_id_p; row_id_p = nullptr;

  NatBC_G( curr_time, dt );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_G(ii);

  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_Solid_FEM::Assem_Tangent_Residual(
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres )
{
  const int nElem = locelem->get_nlocalele();

  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_dot_v = dot_velo -> GetLocalArray();
  const std::vector<double> array_dot_p = dot_pres -> GetLocalArray();

  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  PetscInt * row_id_u = new PetscInt [3*nLocBas];
  PetscInt * row_id_p = new PetscInt [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    const std::vector<int> IEN = locien -> get_LIEN( ee );
    fnode -> get_ctrlPts_xyz(nLocBas, &IEN[0], ectrl_x, ectrl_y, ectrl_z);

    const std::vector<double> local_dot_d = GetLocal( array_dot_d, IEN, nLocBas, 3 );
    const std::vector<double> local_dot_v = GetLocal( array_dot_v, IEN, nLocBas, 3 );
    const std::vector<double> local_dot_p = GetLocal( array_dot_p, IEN, nLocBas, 1 );

    const std::vector<double> local_d = GetLocal( array_d, IEN, nLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, IEN, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN, nLocBas, 1 );

    locassem -> Assem_Tangent_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0],
        &local_dot_p[0], &local_d[0], &local_v[0], &local_p[0],
        ectrl_x, ectrl_y, ectrl_z, &zero_prestress[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = IEN[ii];
      row_id_p[ii] = dof_mat * nbc -> get_LID(0, loc_index);

      for(int mm=0; mm<3; ++mm)
        row_id_u[3*ii + mm] = dof_mat * nbc -> get_LID(mm+1, loc_index) + mm + 1;
    }

    MatSetValues(K, 3*nLocBas, row_id_u, 3*nLocBas, row_id_u,
        locassem->Tangent00, ADD_VALUES);

    MatSetValues(K, 3*nLocBas, row_id_u,   nLocBas, row_id_p,
        locassem->Tangent01, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_u,
        locassem->Tangent10, ADD_VALUES);

    MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p,
        locassem->Tangent11, ADD_VALUES);

    VecSetValues(G, 3*nLocBas, row_id_u, locassem->Residual0, ADD_VALUES);
    VecSetValues(G,   nLocBas, row_id_p, locassem->Residual1, ADD_VALUES);
  }

  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_id_u; row_id_u = nullptr;
  delete [] row_id_p; row_id_p = nullptr;

  NatBC_G( curr_time, dt );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG(ii);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

// EOF
