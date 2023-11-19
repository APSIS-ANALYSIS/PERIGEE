#include "PGAssem_FSI.hpp"

PGAssem_FSI::PGAssem_FSI( 
    IPLocAssem_2x2Block * const &locassem_f_ptr,
    IPLocAssem_2x2Block * const &locassem_s_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_v,
    const ALocal_IEN * const &aien_p,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const ALocal_NBC * const &part_nbc_v,
    const ALocal_NBC * const &part_nbc_p,
    const ALocal_EBC * const &part_ebc,
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( locassem_f_ptr->get_nLocBas_0() ), 
  snLocBas( locassem_f_ptr->get_snLocBas_0() ),
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
  
  // Create matrix with routh preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 4*in_nz_estimate, PETSC_NULL, 4*in_nz_estimate, PETSC_NULL, &K);

  // Create vector
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();
  
  Assem_nonzero_estimate( alelem_ptr, locassem_f_ptr, locassem_s_ptr, elements, quads,
     aien_v, aien_p, pnode_v, part_nbc_v, part_nbc_p, part_ebc, gbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  // Get the maximum nonzeros in each processor
  int Kdnz_max = *std::max_element( Kdnz.begin(), Kdnz.end() );
  int Konz_max = *std::max_element( Konz.begin(), Konz.end() );

  // Reduce the values to processor 0 and print on screen
  MPI_Reduce( SYS_T::get_MPI_rank() ? &Kdnz_max : MPI_IN_PLACE, &Kdnz_max, 1, MPI_INT,
      MPI_MAX, 0, PETSC_COMM_WORLD );

  MPI_Reduce( SYS_T::get_MPI_rank() ? &Konz_max : MPI_IN_PLACE, &Konz_max, 1, MPI_INT,
      MPI_MAX, 0, PETSC_COMM_WORLD );

  if( SYS_T::get_MPI_rank() == 0 ) 
    std::cout<<"===> K dnz max: "<<Kdnz_max<<", and K onz max: "<<Konz_max<<std::endl;

  MPI_Barrier(PETSC_COMM_WORLD);

  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_FSI::~PGAssem_FSI()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_FSI::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const APart_Node * const &pnode_v,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const int nElem = alelem_ptr->get_nlocalele();

  lassem_f_ptr->Assem_Estimate();
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

    if( alelem_ptr->get_elem_tag(ee) == 0 )
    {
      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_f_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K, 3*nLocBas, row_id_v,   nLocBas, row_id_p, lassem_f_ptr->Tangent01, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_v, lassem_f_ptr->Tangent10, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_f_ptr->Tangent11, ADD_VALUES);
    }
    else
    {
      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K, 3*nLocBas, row_id_v,   nLocBas, row_id_p, lassem_s_ptr->Tangent01, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent10, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_s_ptr->Tangent11, ADD_VALUES);
    }
  }

  delete [] row_id_v; row_id_v = nullptr; delete [] row_id_p; row_id_p = nullptr;

  // Resis BC for K and G
  PDNSolution * temp = new PDNSolution_V( pnode_v, 0, false, "aux" );

  NatBC_Resis_KG( 0.1, 0.1, temp, temp, temp, lassem_f_ptr, elements, quad_s,nbc_v, ebc_part, gbc );

  delete temp; temp = nullptr;

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_KG( nbc_v, nbc_p );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_FSI::Assem_mass_residual(
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const Tissue_prestress * const &ps_ptr, 
    const Tissue_property * const &tp_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();

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

    for(int ii=0; ii<nLocBas; ++ii)
    {  
      row_id_v[3*ii  ] = nbc_v -> get_LID(0, IEN_v[ii]);
      row_id_v[3*ii+1] = nbc_v -> get_LID(1, IEN_v[ii]);
      row_id_v[3*ii+2] = nbc_v -> get_LID(2, IEN_v[ii]);
    }

    for(int ii=0; ii<nLocBas; ++ii) 
      row_id_p[ii] = nbc_p -> get_LID( IEN_p[ii] );

    const std::vector<double> local_d = GetLocal( array_d, IEN_v, nLocBas, 3 ); 
    const std::vector<double> local_v = GetLocal( array_v, IEN_v, nLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, IEN_p, nLocBas, 1 );

    if( alelem_ptr->get_elem_tag(ee) == 0 )
    {
      lassem_f_ptr->Assem_Mass_Residual(&local_d[0], &local_v[0], &local_p[0], elementv, 
          ectrl_x, ectrl_y, ectrl_z, quad_v);

      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_f_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_f_ptr->Tangent11, ADD_VALUES);

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_f_ptr->Residual0, ADD_VALUES);
    }
    else
    {
      // For solid element, quaprestress will return a vector of length nqp x 6
      // for the prestress values at the quadrature points
      const std::vector<double> quaprestress = ps_ptr->get_prestress( ee );

      // For solid element, the direction basis includes 3 Vector_3 for each node.
      // ebasis_r, ebasis_l, and ebasis_c are vectors of length nLocBas.
      std::vector<Vector_3> ebasis_r(nLocBas);
      std::vector<Vector_3> ebasis_l(nLocBas);
      std::vector<Vector_3> ebasis_c(nLocBas);

      for(int ii=0; ii<nLocBas; ++ii)
      {
	ebasis_r[ii] = tp_ptr -> get_basis_r(IEN_v[ii]);
	ebasis_l[ii] = tp_ptr -> get_basis_l(IEN_v[ii]);
	ebasis_c[ii] = tp_ptr -> get_basis_c(IEN_v[ii]);
      }

      lassem_s_ptr->Assem_Mass_Residual(&local_d[0], &local_v[0], &local_p[0], elementv, 
          ectrl_x, ectrl_y, ectrl_z, &quaprestress[0], quad_v, ebasis_r, ebasis_l, ebasis_c);

      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_s_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_s_ptr->Tangent11, ADD_VALUES);

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_s_ptr->Residual0, ADD_VALUES);
    } 

  }

  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  ectrl_x = nullptr; ectrl_y = nullptr; ectrl_z = nullptr;
  delete [] row_id_v; delete [] row_id_p;
  row_id_v = nullptr; row_id_p = nullptr;

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_KG( nbc_v, nbc_p );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_FSI::Assem_Residual(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    const PDNSolution * const &dot_velo_np1,
    const PDNSolution * const &velo_np1,
    const PDNSolution * const &disp_np1,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const Tissue_prestress * const &ps_ptr )
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

    if( alelem_ptr->get_elem_tag(ee) == 0 )
    {
      lassem_f_ptr -> Assem_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0], &local_dot_p[0],
          &local_d[0], &local_v[0], &local_p[0], elementv, ectrl_x, ectrl_y, ectrl_z, quad_v );

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_f_ptr->Residual0, ADD_VALUES);
      VecSetValues(G,   nLocBas, row_id_p, lassem_f_ptr->Residual1, ADD_VALUES);
    }
    else
    {
      // For solid element, quaprestress will return a vector of length nqp x 6
      // for the prestress values at the quadrature points
      const std::vector<double> quaprestress = ps_ptr->get_prestress( ee );

      lassem_s_ptr -> Assem_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0], &local_dot_p[0],
          &local_d[0], &local_v[0], &local_p[0], elementv, ectrl_x, ectrl_y, ectrl_z, &quaprestress[0], quad_v );

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_s_ptr->Residual0, ADD_VALUES);
      VecSetValues(G,   nLocBas, row_id_p, lassem_s_ptr->Residual1, ADD_VALUES);
    }
  }

  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  ectrl_x = nullptr; ectrl_y = nullptr; ectrl_z = nullptr;
  delete [] row_id_v; delete [] row_id_p;
  row_id_v = nullptr; row_id_p = nullptr;

  // Backflow stabilization
  BackFlow_G( dot_disp, disp, velo, lassem_f_ptr, elements, quad_s, nbc_v, ebc_part );

  // Resistance BC for G
  NatBC_Resis_G( curr_time, dt, disp_np1, dot_velo_np1, velo_np1, lassem_f_ptr, elements, quad_s,
      nbc_v, ebc_part, gbc );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_G( nbc_v, nbc_p );

  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_FSI::Assem_Tangent_Residual(
    const double &curr_time, const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &dot_pres,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    const PDNSolution * const &dot_velo_np1,
    const PDNSolution * const &velo_np1,
    const PDNSolution * const &disp_np1,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    IPLocAssem_2x2Block * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &fnode_ptr,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const Tissue_prestress * const &ps_ptr )
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

    if( alelem_ptr->get_elem_tag(ee) == 0 )
    {
      lassem_f_ptr -> Assem_Tangent_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0], &local_dot_p[0],
          &local_d[0], &local_v[0], &local_p[0], elementv, ectrl_x, ectrl_y, ectrl_z, quad_v );

      MatSetValues(K, 3*nLocBas, row_id_v, 3*nLocBas, row_id_v, lassem_f_ptr->Tangent00, ADD_VALUES);

      MatSetValues(K, 3*nLocBas, row_id_v,   nLocBas, row_id_p, lassem_f_ptr->Tangent01, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p, 3*nLocBas, row_id_v, lassem_f_ptr->Tangent10, ADD_VALUES);

      MatSetValues(K,   nLocBas, row_id_p,   nLocBas, row_id_p, lassem_f_ptr->Tangent11, ADD_VALUES);

      VecSetValues(G, 3*nLocBas, row_id_v, lassem_f_ptr->Residual0, ADD_VALUES);
      VecSetValues(G,   nLocBas, row_id_p, lassem_f_ptr->Residual1, ADD_VALUES);
    }
    else
    {
      // For solid element, quaprestress will return a vector of length nqp x 6
      // for the prestress values at the quadrature points
      const std::vector<double> quaprestress = ps_ptr->get_prestress( ee );

      lassem_s_ptr -> Assem_Tangent_Residual( curr_time, dt, &local_dot_d[0], &local_dot_v[0], &local_dot_p[0],
          &local_d[0], &local_v[0], &local_p[0], elementv, ectrl_x, ectrl_y, ectrl_z, &quaprestress[0], quad_v );

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

  // Backflow stabilization
  BackFlow_KG( dt, dot_disp, disp, velo, lassem_f_ptr, elements, quad_s, nbc_v, ebc_part );

  // Resistance BC for G
  NatBC_Resis_KG( curr_time, dt, disp_np1, dot_velo_np1, velo_np1, lassem_f_ptr, elements, quad_s,
      nbc_v, ebc_part, gbc );

  VecAssemblyBegin(G); VecAssemblyEnd(G);

  EssBC_KG( nbc_v, nbc_p );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G); VecAssemblyEnd(G);
}

void PGAssem_FSI::EssBC_KG( const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p )
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


void PGAssem_FSI::EssBC_G( const ALocal_NBC * const &nbc_v, 
    const ALocal_NBC * const &nbc_p )
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

double PGAssem_FSI::Assem_surface_flowrate(
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    const std::vector<int> LSIEN = ebc_part -> get_SIEN( ebc_id, ee );

    ebc_part -> get_ctrlPts_xyz( ebc_id, ee, sctrl_x, sctrl_y, sctrl_z );
    
    const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, LSIEN, snLocBas, 3 );
   
    esum += lassem_ptr -> get_flowrate( &local_d[0], &local_v[0], element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum = 0.0;

  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}

double PGAssem_FSI::Assem_surface_flowrate(
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_InflowBC * const &infbc_part,
    const int &nbc_id )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = infbc_part -> get_num_local_cell(nbc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    const std::vector<int> LSIEN = infbc_part -> get_SIEN( nbc_id, ee );
    
    infbc_part -> get_ctrlPts_xyz( nbc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 );
    const std::vector<double> local_v = GetLocal( array_v, LSIEN, snLocBas, 3 );

    esum += lassem_ptr -> get_flowrate( &local_d[0], &local_v[0], element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum = 0.0;

  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}

double PGAssem_FSI::Assem_surface_ave_pressure(
    const PDNSolution * const &disp,
    const PDNSolution * const &pres,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_v,
    const ALocal_EBC * const &ebc_p,
    const int &ebc_id )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = ebc_v -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {

    const std::vector<int> LSIEN_v = ebc_v -> get_SIEN( ebc_id, ee );
    const std::vector<int> LSIEN_p = ebc_p -> get_SIEN( ebc_id, ee );

    ebc_v -> get_ctrlPts_xyz( ebc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    const std::vector<double> local_d = GetLocal( array_d, LSIEN_v, snLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, LSIEN_p, snLocBas, 1 );

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( &local_d[0], &local_p[0], element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);
    
    val_pres += ele_pres;
    val_area += ele_area; 
  }

  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}

double PGAssem_FSI::Assem_surface_ave_pressure(
    const PDNSolution * const &disp,
    const PDNSolution * const &pres,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_InflowBC * const &infbc_part,
    const int &nbc_id )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_p = pres -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = infbc_part -> get_num_local_cell(nbc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    const std::vector<int> LSIEN = infbc_part -> get_SIEN( nbc_id, ee );

    infbc_part -> get_ctrlPts_xyz( nbc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 );
    const std::vector<double> local_p = GetLocal( array_p, LSIEN, snLocBas, 1 );

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( &local_d[0], &local_p[0], element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}

void PGAssem_FSI::NatBC_G( const double &curr_time, const double &dt,
    const PDNSolution * const &disp,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part )
{
  const std::vector<double> array_d = disp -> GetLocalArray(); 

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscInt * srow_index = new PetscInt [3*snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      const std::vector<int> LSIEN = ebc_part -> get_SIEN(ebc_id, ee); 
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 ); 

      lassem_f_ptr -> Assem_Residual_EBC( ebc_id, curr_time, dt, &local_d[0], element_s,
          sctrl_x, sctrl_y, sctrl_z, quad_s );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_index[3*ii+0] = nbc_v -> get_LID(0, LSIEN[ii]);
        srow_index[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
        srow_index[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
      }

      VecSetValues(G, 3*snLocBas, srow_index, lassem_f_ptr->sur_Residual0, ADD_VALUES); 
    }
  }

  delete [] srow_index; srow_index = nullptr;
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

void PGAssem_FSI::NatBC_Resis_G( const double &curr_time, const double &dt,
    const PDNSolution * const &disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const std::vector<double> array_d = disp -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscInt * srow_index = new PetscInt [3*snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id
    const double dot_flrate = Assem_surface_flowrate( disp, dot_velo, lassem_f_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Calculate flow rate for face with ebc_id
    const double flrate = Assem_surface_flowrate( disp, velo, lassem_f_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the pressure value on the outlet surfaces
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate, curr_time + dt );

    // P_n+alpha_f
    const double val = P_n + lassem_f_ptr->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      const std::vector<int> LSIEN = ebc_part -> get_SIEN( ebc_id, ee );

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 );

      lassem_f_ptr->Assem_Residual_EBC_Resistance( val, &local_d[0],
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_index[3*ii+0] = nbc_v -> get_LID(0, LSIEN[ii]);
        srow_index[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
        srow_index[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
      }

      VecSetValues(G, 3*snLocBas, srow_index, lassem_f_ptr->sur_Residual0, ADD_VALUES); 
    }
  }

  delete [] srow_index; srow_index = nullptr;
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

void PGAssem_FSI::NatBC_Resis_KG( const double &curr_time, const double &dt,
    const PDNSolution * const &disp,
    const PDNSolution * const &dot_velo,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  PetscScalar * Tan;
  PetscInt * scol_idx;
  Vector_3 out_n;
  std::vector<double> intNB;
  std::vector<int> map_Bj;

  const double a_f = lassem_f_ptr->get_model_para_1();

  const double a_gamma = lassem_f_ptr->get_model_para_2();

  const double dd_dv = dt * a_f * a_gamma;

  const std::vector<double> array_d = disp -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscScalar * Res = new PetscScalar [3*snLocBas];
  PetscInt * srow_idx = new PetscInt [3*snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id
    const double dot_flrate = Assem_surface_flowrate( disp, dot_velo, lassem_f_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Calculate flow rate for face with ebc_id
    const double flrate = Assem_surface_flowrate( disp, velo, lassem_f_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the pressure value on the outlet surfaces
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate, curr_time + dt );

    // P_n+alpha_f
    const double resis_val = P_n + a_f * (P_np1 - P_n);

    // Get m := dP/dQ
    const double m_val = gbc -> get_m( ebc_id, dot_flrate, flrate );

    // Get n := dP/d(dot_Q)
    const double n_val = gbc -> get_n( ebc_id, dot_flrate, flrate );

    // Define alpha_f x n + alpha_f x gamma x dt x m
    const double coef = a_f * n_val + dd_dv * m_val;

    const int num_face_nodes = ebc_part -> get_num_face_nodes(ebc_id);

    if(num_face_nodes > 0)
    {
      Tan = new PetscScalar [snLocBas * 3 * num_face_nodes * 3];
      scol_idx = new PetscInt [num_face_nodes * 3];
      out_n = ebc_part -> get_outvec( ebc_id );
      intNB = ebc_part -> get_intNA( ebc_id );
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
      const std::vector<int> LSIEN = ebc_part -> get_SIEN( ebc_id, ee );

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      const std::vector<double> local_d = GetLocal( array_d, LSIEN, snLocBas, 3 );

      lassem_f_ptr->Assem_Residual_EBC_Resistance( 1.0, &local_d[0],
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii  ] = resis_val * lassem_f_ptr->sur_Residual0[3*ii];
        Res[3*ii+1] = resis_val * lassem_f_ptr->sur_Residual0[3*ii+1];
        Res[3*ii+2] = resis_val * lassem_f_ptr->sur_Residual0[3*ii+2];

        srow_idx[3*ii  ] = nbc_v -> get_LID(0, LSIEN[ii]);
        srow_idx[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
        srow_idx[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            Tan[temp_row + 3*B + 0] = coef * lassem_f_ptr->sur_Residual0[3*A+ii] * intNB[B] * out_n.x();
            Tan[temp_row + 3*B + 1] = coef * lassem_f_ptr->sur_Residual0[3*A+ii] * intNB[B] * out_n.y();
            Tan[temp_row + 3*B + 2] = coef * lassem_f_ptr->sur_Residual0[3*A+ii] * intNB[B] * out_n.z();
          }
        }
      }

      for(int ii=0; ii<num_face_nodes; ++ii)
      {
        scol_idx[ii*3+0] = map_Bj[ii*3+0];
        scol_idx[ii*3+1] = map_Bj[ii*3+1];
        scol_idx[ii*3+2] = map_Bj[ii*3+2];
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
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

void PGAssem_FSI::BackFlow_G( 
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscInt * srow_index = new PetscInt [3*snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      const std::vector<int> LSIEN = ebc_part -> get_SIEN(ebc_id, ee);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      const std::vector<double> local_d     = GetLocal( array_d,     LSIEN, snLocBas, 3 );
      const std::vector<double> local_dot_d = GetLocal( array_dot_d, LSIEN, snLocBas, 3 );
      const std::vector<double> local_v     = GetLocal( array_v,     LSIEN, snLocBas, 3 );  

      lassem_f_ptr->Assem_Residual_BackFlowStab( &local_dot_d[0], &local_d[0], &local_v[0],
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_index[3*ii+0] = nbc_v -> get_LID(0, LSIEN[ii]);
        srow_index[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
        srow_index[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
      }

      VecSetValues(G, 3*snLocBas, srow_index, lassem_f_ptr->sur_Residual0, ADD_VALUES); 
    }
  }

  delete [] srow_index; srow_index = nullptr;
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

void PGAssem_FSI::BackFlow_KG( const double &dt,
    const PDNSolution * const &dot_disp,
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NBC * const &nbc_v,
    const ALocal_EBC * const &ebc_part )
{
  const std::vector<double> array_d = disp -> GetLocalArray();
  const std::vector<double> array_dot_d = dot_disp -> GetLocalArray();
  const std::vector<double> array_v = velo -> GetLocalArray();

  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  PetscInt * srow_index = new PetscInt [3*snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      const std::vector<int> LSIEN = ebc_part -> get_SIEN(ebc_id, ee);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      const std::vector<double> local_d     = GetLocal( array_d,     LSIEN, snLocBas, 3 );
      const std::vector<double> local_dot_d = GetLocal( array_dot_d, LSIEN, snLocBas, 3 );
      const std::vector<double> local_v     = GetLocal( array_v,     LSIEN, snLocBas, 3 );

      lassem_f_ptr->Assem_Tangent_Residual_BackFlowStab( dt, &local_dot_d[0], &local_d[0], 
          &local_v[0], element_s, sctrl_x, sctrl_y, sctrl_z, quad_s );

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_index[3*ii+0] = nbc_v -> get_LID(0, LSIEN[ii]);
        srow_index[3*ii+1] = nbc_v -> get_LID(1, LSIEN[ii]);
        srow_index[3*ii+2] = nbc_v -> get_LID(2, LSIEN[ii]);
      }

      VecSetValues(G, 3*snLocBas, srow_index, lassem_f_ptr->sur_Residual0, ADD_VALUES);

      MatSetValues(K, 3*snLocBas, srow_index, 3*snLocBas, srow_index,
          lassem_f_ptr->sur_Tangent00, ADD_VALUES);
    }
  }

  delete [] srow_index; srow_index = nullptr;
  delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;
}

// EOF
