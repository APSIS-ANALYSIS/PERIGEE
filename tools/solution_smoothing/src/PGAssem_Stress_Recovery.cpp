#include "PGAssem_Stress_Recovery.hpp"

PGAssem_Stress_Recovery::PGAssem_Stress_Recovery(
    IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ), nlgn( pnode_ptr->get_nlocghonode() )
{
  const int nlocrow = 6 * pnode_ptr -> get_nlocalnode();

  // Allocate the sparse matrix K
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 3*in_nz_estimate, NULL, 3*in_nz_estimate, NULL, &K);

  // Allocate the vector G
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, aien_ptr, pnode_ptr );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);
  
  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_Stress_Recovery::~PGAssem_Stress_Recovery()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_Stress_Recovery::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &pnode_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();
  
  lassem_ptr->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * 6];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index = lien_ptr -> get_LIEN(ee, ii);
      
      for(int mm=0; mm<6; ++mm)
        row_index[6 * ii + mm] = 6 * pnode_ptr->get_local_to_global(loc_index) + mm;
    }

    MatSetValues(K, 6 * nLocBas, row_index, 6 * nLocBas, row_index, lassem_ptr->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Stress_Recovery::Assem_residual(
    const PDNSolution * const &isol,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    const IQuadPts * const &quad_v,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const APart_Node * const &pnode_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();

  double * array_a = new double [nlgn * 3];
  double * local_a = new double [nLocBas * 3];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * 6];

  isol -> GetLocalArray( array_a );
  
  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr -> get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Residual(local_a, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<6; ++mm)
        row_index[6*ii+mm] = 6 * pnode_ptr->get_local_to_global(IEN_e[ii]) + mm;
    }

    VecSetValues(G, 6 * nLocBas, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Stress_Recovery::Assem_mass_residual(
    const PDNSolution * const &isol,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    const IQuadPts * const &quad_v,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &fnode_ptr,
    const APart_Node * const &pnode_ptr )
{
  const int nElem = alelem_ptr->get_nlocalele();

  double * array_a = new double [nlgn * 3];
  double * local_a = new double [nLocBas * 3];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * 6];
  isol->GetLocalArray( array_a );
  
  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    
    GetLocal(array_a, IEN_e, local_a);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Mass_Residual( local_a, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<6; ++mm)
        row_index[6*ii+mm] = 6 * pnode_ptr->get_local_to_global(IEN_e[ii]) + mm;
    }

    MatSetValues(K, 6 * nLocBas, row_index, 6 * nLocBas, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, 6 * nLocBas, row_index, lassem_ptr->Residual, ADD_VALUES);
  }
  
  delete [] array_a; array_a = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

//EOF