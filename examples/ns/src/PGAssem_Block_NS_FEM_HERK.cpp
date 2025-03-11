#include "PGAssem_Block_NS_FEM_HERK.hpp"

PGAssem_Block_NS_FEM_HERK::PGAssem_Block_NS_FEM_HERK(
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    std::unique_ptr<ALocal_NBC> in_nbc,
    std::unique_ptr<ALocal_EBC> in_ebc,
    std::unique_ptr<PLocAssem_Block_VMS_NS_HERK> in_locassem,    
    const int &in_nz_estimate )
: locien( std::move(in_locien) ),
  locelem( std::move(in_locelem) ),
  fnode( std::move(in_fnode) ),
  pnode( std::move(in_pnode) ),
  nbc( std::move(in_nbc) ),
  ebc( std::move(in_ebc) ),
  locassem(std::move(in_locassem)),
  nLocBas( locassem->get_nLocBas() ),
  snLocBas( locassem->get_snLocBas() ),
  dof_sol( pnode->get_dof() ),
  dof_mat_v( locassem->get_dof_mat_0() ),
  dof_mat_p( locassem->get_dof_mat_1() ),
  num_ebc( ebc->get_num_ebc() ),
  nlgn( pnode->get_nlocghonode() )
{
  SYS_T::print_fatal_if(dof_sol != locassem->get_dof(),
      "PGAssem_Block_NS_FEM_HERK::dof_sol != locassem->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat_v + dof_mat_p != nbc->get_dof_LID(),
      "PGAssem_Block_NS_FEM_HERK::dof_mat != nbc->get_dof_LID(). \n");
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_Block_NS_FEM_HERK, snLocBas has to be uniform. \n");
  }

  const int nlocrow_v  = dof_mat_v * pnode->get_nlocalnode(); 
  const int nlocrow_p  = dof_mat_p * pnode->get_nlocalnode();

  // Allocate the block matrix K
  // D matrix
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL,
      dof_mat_p*in_nz_estimate, NULL, &subK[0]);

  // C matrix
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_v, PETSC_DETERMINE,
    PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL, 
    dof_mat_v*in_nz_estimate, NULL, &subK[1]);

  // B matrix
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_p, PETSC_DETERMINE,
    PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
    dof_mat_p*in_nz_estimate, NULL, &subK[2]);

  // A matrix
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
      dof_mat_v*in_nz_estimate, NULL, &subK[3]);

  // A_tilde matrix
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_v, PETSC_DETERMINE,
    PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
    dof_mat_v*in_nz_estimate, NULL, &subK[4]);

  // Allocate the sub-vectors
  VecCreate(PETSC_COMM_WORLD, &subG[0]);
  VecCreate(PETSC_COMM_WORLD, &subG[1]);

  VecSetSizes(subG[0], nlocrow_p, PETSC_DECIDE);
  VecSetSizes(subG[1], nlocrow_v, PETSC_DECIDE);

  VecSetFromOptions(subG[0]);
  VecSetFromOptions(subG[1]);

  VecSet(subG[0], 0.0);
  VecSet(subG[1], 0.0);
  
  VecSetOption(subG[0], VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(subG[1], VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate();

  // Reallocate matrices with precise preallocation
  std::vector<int> Kdnz, Konz;

  PETSc_T::Get_dnz_onz(subK[0], Kdnz, Konz);
  MatDestroy(&subK[0]);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &subK[0]);
  
  PETSc_T::Get_dnz_onz(subK[1], Kdnz, Konz);
  MatDestroy(&subK[1]);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &subK[1]);
  
  PETSc_T::Get_dnz_onz(subK[2], Kdnz, Konz);
  MatDestroy(&subK[2]);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &subK[2]);

  PETSc_T::Get_dnz_onz(subK[3], Kdnz, Konz);
  MatDestroy(&subK[3]);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &subK[3]);
 
  PETSc_T::Get_dnz_onz(subK[4], Kdnz, Konz);
  MatDestroy(&subK[4]);
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &subK[4]);
}

PGAssem_Block_NS_FEM_HERK::~PGAssem_Block_NS_FEM_HERK()
{
  VecDestroy(&subG[0]); VecDestroy(&subG[1]); 
  MatDestroy(&subK[0]); MatDestroy(&subK[1]);
  MatDestroy(&subK[2]); MatDestroy(&subK[3]);
  MatDestroy(&subK[4]);
}

void PGAssem_Block_NS_FEM_HERK::EssBC_KG()
{
  // pressure dof comes from field 0, to be inserted in subK[0] and subG[0]
  const int local_dir = nbc->get_Num_LD(0);

  if( local_dir > 0 )
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc->get_LDN(0, i) * dof_mat_p;
      
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
      MatSetValue(subK[0], row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(0);
  if( local_sla > 0 )
  {
    for(int i=0; i<local_sla; ++i)
    {
      const int row = nbc->get_LPSN(0, i) * dof_mat_p;
      const int col = nbc->get_LPMN(0, i) * dof_mat_p;
      MatSetValue(subK[0], row, col, 1.0, ADD_VALUES);
      MatSetValue(subK[0], row, row, -1.0, ADD_VALUES);
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  // velocity dofs
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc->get_Num_LD(field);

    if( local_dir > 0 )
    {
      for(int i=0; i<local_dir; ++i)
      {
        const int row = nbc->get_LDN(field, i) * dof_mat_v + field - 1;
        MatSetValue(subK[3], row, row, 1.0, ADD_VALUES);
        MatSetValue(subK[4], row, row, 0.0, ADD_VALUES);
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);        
      }
    }

    const int local_sla = nbc->get_Num_LPS(field);

    if( local_sla > 0 )
    {
      for(int i=0; i<local_sla; ++i)
      {
        const int row = nbc->get_LPSN(field, i) * dof_mat_v + field - 1;
        const int col = nbc->get_LPMN(field, i) * dof_mat_v + field - 1;
        MatSetValue(subK[3], row, col, 1.0, ADD_VALUES);
        MatSetValue(subK[3], row, row, -1.0, ADD_VALUES);
        MatSetValue(subK[4], row, row, 0.0, ADD_VALUES);
        MatSetValue(subK[4], row, col, 0.0, ADD_VALUES);
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);        
      }
    }
  }
}

void PGAssem_Block_NS_FEM_HERK::EssBC_G()
{
  // pres field is 0, to be inserted to subG[0]
  const int local_dir = nbc->get_Num_LD(0);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc->get_LDN(0, ii) * dof_mat_p;
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(0);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc->get_LPSN(0, ii) * dof_mat_p;
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  // velo fields from 1 to 3
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc->get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc->get_LDN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc->get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc->get_LPSN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
      }
    }
  }
}

void PGAssem_Block_NS_FEM_HERK::Assem_nonzero_estimate()
{
  const int nElem = locelem->get_nlocalele();
  const int loc_dof_v = dof_mat_v * nLocBas;
  const int loc_dof_p = dof_mat_p * nLocBas;
  
  locassem->Assem_Estimate();

  PetscInt * row_idx_v = new PetscInt [nLocBas * dof_mat_v];
  PetscInt * row_idx_p = new PetscInt [nLocBas * dof_mat_p];

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int loc_index  = locien->get_LIEN(ee, ii);

      row_idx_v[3*ii]   = 3 * nbc->get_LID( 1, loc_index );
      row_idx_v[3*ii+1] = 3 * nbc->get_LID( 2, loc_index ) + 1;
      row_idx_v[3*ii+2] = 3 * nbc->get_LID( 3, loc_index ) + 2;
      
      row_idx_p[ii] = nbc->get_LID( 0, loc_index );
    }

    MatSetValues(subK[0], loc_dof_p, row_idx_p, loc_dof_p, row_idx_p, locassem->Tangent0, ADD_VALUES);
    MatSetValues(subK[1], loc_dof_p, row_idx_p, loc_dof_v, row_idx_v, locassem->Tangent1, ADD_VALUES);
    MatSetValues(subK[2], loc_dof_v, row_idx_v, loc_dof_p, row_idx_p, locassem->Tangent2, ADD_VALUES);
    MatSetValues(subK[3], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent3, ADD_VALUES);    
    MatSetValues(subK[4], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent4, ADD_VALUES); 
  }

  delete [] row_idx_v; row_idx_v = nullptr;
  delete [] row_idx_p; row_idx_p = nullptr;

  // Create a temporary zero solution vector to feed Natbc_Resis_KG
  // PDNSolution * temp = new PDNSolution_NS( pnode.get(), 0, false );

  // // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  // NatBC_Resis_KG( 0.0, 0.1, temp, temp );

  // delete temp;

  VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
  VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);

  EssBC_KG();

  MatAssemblyBegin(subK[0], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[0], MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(subK[1], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[1], MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(subK[2], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[2], MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(subK[3], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[3], MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(subK[4], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[4], MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
  VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);
}

void PGAssem_Block_NS_FEM_HERK::Assem_tangent_residual_substep(
  const int &substep_index,
  PDNSolution ** const &cur_velo_sols,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution ** const &pre_velo_sols,
  PDNSolution * const &pre_velo,
  PDNSolution ** const &pre_pres_sols,
  PDNSolution * const &pre_velo_before,   
  const ITimeMethod_RungeKutta * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = locelem->get_nlocalele();
const int loc_dof_v = dof_mat_v * nLocBas;
const int loc_dof_p = dof_mat_p * nLocBas;

std::vector<std::vector<double>> array_cur_velo_sols(substep_index+1);
std::vector<std::vector<double>> array_cur_pres_sols(substep_index+1);
std::vector<std::vector<double>> array_pre_velo_sols(tm_RK_ptr->get_RK_step());
std::vector<std::vector<double>> array_pre_pres_sols(tm_RK_ptr->get_RK_step());

for(int ii = 0; ii < substep_index+1; ++ii)
{
  array_cur_velo_sols[ii] = cur_velo_sols[ii] -> GetLocalArray();
  array_cur_pres_sols[ii] = cur_pres_sols[ii] -> GetLocalArray();
}

for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
{
  array_pre_velo_sols[ii] = pre_velo_sols[ii] -> GetLocalArray();
  array_pre_pres_sols[ii] = pre_pres_sols[ii] -> GetLocalArray();
}

const std::vector<double> array_pre_velo = pre_velo -> GetLocalArray();
const std::vector<double> array_pre_velo_before = pre_velo_before -> GetLocalArray();

double * ectrl_x = new double [nLocBas];
double * ectrl_y = new double [nLocBas];
double * ectrl_z = new double [nLocBas];
PetscInt * row_idx_v = new PetscInt [nLocBas * dof_mat_v];
PetscInt * row_idx_p = new PetscInt [nLocBas * dof_mat_p];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = locien->get_LIEN(ee);

  std::vector<std::vector<double>> local_cur_velo_sols(substep_index+1);    
  std::vector<std::vector<double>> local_cur_pres_sols(substep_index+1); 
  std::vector<std::vector<double>> local_pre_velo_sols(tm_RK_ptr->get_RK_step()); 
  std::vector<std::vector<double>> local_pre_pres_sols(tm_RK_ptr->get_RK_step()); 

  for(int ii = 0; ii < substep_index+1; ++ii)
  {
    local_cur_velo_sols[ii] = GetLocal( array_cur_velo_sols[ii], IEN_e, nLocBas, 3 );
    local_cur_pres_sols[ii] = GetLocal( array_cur_pres_sols[ii], IEN_e, nLocBas, 1 );
  }

  for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
  {
    local_pre_velo_sols[ii] = GetLocal( array_pre_velo_sols[ii], IEN_e, nLocBas, 3 );
    local_pre_pres_sols[ii] = GetLocal( array_pre_pres_sols[ii], IEN_e, nLocBas, 1 );
  }

  const std::vector<double> local_pre_velo = GetLocal( array_pre_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_pre_velo_before = GetLocal( array_pre_velo_before, IEN_e, nLocBas, 3 );

  fnode->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  locassem->Assem_Tangent_Residual_Sub(curr_time, dt, substep_index, tm_RK_ptr, local_cur_velo_sols, local_cur_pres_sols,
      local_pre_velo_sols, local_pre_pres_sols, local_pre_velo, local_pre_velo_before, ectrl_x, ectrl_y, ectrl_z);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    const int loc_index  = locien->get_LIEN(ee, ii);

    row_idx_v[3*ii]   = 3 * nbc->get_LID( 1, loc_index );
    row_idx_v[3*ii+1] = 3 * nbc->get_LID( 2, loc_index ) + 1;
    row_idx_v[3*ii+2] = 3 * nbc->get_LID( 3, loc_index ) + 2;

    row_idx_p[ii] = nbc->get_LID( 0, loc_index );
  }

  MatSetValues(subK[0], loc_dof_p, row_idx_p, loc_dof_p, row_idx_p, locassem->Tangent0, ADD_VALUES);
  MatSetValues(subK[1], loc_dof_p, row_idx_p, loc_dof_v, row_idx_v, locassem->Tangent1, ADD_VALUES);
  MatSetValues(subK[2], loc_dof_v, row_idx_v, loc_dof_p, row_idx_p, locassem->Tangent2, ADD_VALUES);
  MatSetValues(subK[3], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent3, ADD_VALUES);
  MatSetValues(subK[4], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent4, ADD_VALUES);

  VecSetValues(subG[0], loc_dof_p, row_idx_p, locassem->Residual0, ADD_VALUES);
  VecSetValues(subG[1], loc_dof_v, row_idx_v, locassem->Residual1, ADD_VALUES);
}

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_idx_v; row_idx_v = nullptr;
delete [] row_idx_p; row_idx_p = nullptr;

NatBC_G_HERK_Sub( curr_time, dt, substep_index, tm_RK_ptr );

VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);

EssBC_KG();

MatAssemblyBegin(subK[0], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[0], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[1], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[1], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[2], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[2], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[3], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[3], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[4], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[4], MAT_FINAL_ASSEMBLY);

VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);
}

void PGAssem_Block_NS_FEM_HERK::Assem_tangent_residual_finalstep(
  PDNSolution ** const &cur_velo_sols,
  PDNSolution * const &cur_velo,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution ** const &pre_velo_sols,
  PDNSolution * const &pre_velo,
  PDNSolution ** const &pre_pres_sols,
  PDNSolution * const &pre_velo_before,    
  const ITimeMethod_RungeKutta * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = locelem->get_nlocalele();
const int loc_dof_v = dof_mat_v * nLocBas;
const int loc_dof_p = dof_mat_p * nLocBas;

std::vector<std::vector<double>> array_cur_velo_sols(tm_RK_ptr->get_RK_step());
std::vector<std::vector<double>> array_cur_pres_sols(tm_RK_ptr->get_RK_step());
std::vector<std::vector<double>> array_pre_velo_sols(tm_RK_ptr->get_RK_step());
std::vector<std::vector<double>> array_pre_pres_sols(tm_RK_ptr->get_RK_step());

for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
{
  array_cur_velo_sols[ii] = cur_velo_sols[ii] -> GetLocalArray();
  array_cur_pres_sols[ii] = cur_pres_sols[ii] -> GetLocalArray();
  array_pre_velo_sols[ii] = pre_velo_sols[ii] -> GetLocalArray();
  array_pre_pres_sols[ii] = pre_pres_sols[ii] -> GetLocalArray();
}

const std::vector<double> array_cur_velo = cur_velo -> GetLocalArray();
const std::vector<double> array_pre_velo = pre_velo -> GetLocalArray();
const std::vector<double> array_pre_velo_before = pre_velo_before -> GetLocalArray();

double * ectrl_x = new double [nLocBas];
double * ectrl_y = new double [nLocBas];
double * ectrl_z = new double [nLocBas];
PetscInt * row_idx_v = new PetscInt [nLocBas * dof_mat_v];
PetscInt * row_idx_p = new PetscInt [nLocBas * dof_mat_p];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = locien->get_LIEN(ee);

  std::vector<std::vector<double>> local_cur_velo_sols(tm_RK_ptr->get_RK_step());    
  std::vector<std::vector<double>> local_cur_pres_sols(tm_RK_ptr->get_RK_step()); 
  std::vector<std::vector<double>> local_pre_velo_sols(tm_RK_ptr->get_RK_step()); 
  std::vector<std::vector<double>> local_pre_pres_sols(tm_RK_ptr->get_RK_step()); 

  for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
  {
    local_cur_velo_sols[ii] = GetLocal( array_cur_velo_sols[ii], IEN_e, nLocBas, 3 );
    local_cur_pres_sols[ii] = GetLocal( array_cur_pres_sols[ii], IEN_e, nLocBas, 1 );
    local_pre_velo_sols[ii] = GetLocal( array_pre_velo_sols[ii], IEN_e, nLocBas, 3 );
    local_pre_pres_sols[ii] = GetLocal( array_pre_pres_sols[ii], IEN_e, nLocBas, 1 );
  }

  const std::vector<double> local_cur_velo = GetLocal( array_cur_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_pre_velo = GetLocal( array_pre_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_pre_velo_before = GetLocal( array_pre_velo_before, IEN_e, nLocBas, 3 );

  fnode->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  locassem->Assem_Tangent_Residual_Final(curr_time, dt, tm_RK_ptr, local_cur_velo_sols, local_cur_velo,
      local_cur_pres_sols, local_pre_velo_sols, local_pre_velo, local_pre_pres_sols, local_pre_velo_before, 
      ectrl_x, ectrl_y, ectrl_z);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    const int loc_index  = locien->get_LIEN(ee, ii);

    row_idx_v[3*ii]   = 3 * nbc->get_LID( 1, loc_index );
    row_idx_v[3*ii+1] = 3 * nbc->get_LID( 2, loc_index ) + 1;
    row_idx_v[3*ii+2] = 3 * nbc->get_LID( 3, loc_index ) + 2;

    row_idx_p[ii] = nbc->get_LID( 0, loc_index );
  }
  
  MatSetValues(subK[0], loc_dof_p, row_idx_p, loc_dof_p, row_idx_p, locassem->Tangent0, ADD_VALUES);
  MatSetValues(subK[1], loc_dof_p, row_idx_p, loc_dof_v, row_idx_v, locassem->Tangent1, ADD_VALUES);
  MatSetValues(subK[2], loc_dof_v, row_idx_v, loc_dof_p, row_idx_p, locassem->Tangent2, ADD_VALUES);
  MatSetValues(subK[3], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent3, ADD_VALUES);
  MatSetValues(subK[4], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent4, ADD_VALUES);

  VecSetValues(subG[0], loc_dof_p, row_idx_p, locassem->Residual0, ADD_VALUES);
  VecSetValues(subG[1], loc_dof_v, row_idx_v, locassem->Residual1, ADD_VALUES);
}

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_idx_v; row_idx_v = nullptr;
delete [] row_idx_p; row_idx_p = nullptr;

NatBC_G_HERK_Final( curr_time, dt, tm_RK_ptr );

VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);

EssBC_KG();

MatAssemblyBegin(subK[0], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[0], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[1], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[1], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[2], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[2], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[3], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[3], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[4], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[4], MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);
}

void PGAssem_Block_NS_FEM_HERK::Assem_tangent_residual_presstage(
  PDNSolution * const &cur_dot_velo,
  PDNSolution ** const &cur_velo_sols,
  PDNSolution * const &cur_velo,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution * const &pre_velo,
  PDNSolution * const &cur_pres,    
  const ITimeMethod_RungeKutta * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = locelem->get_nlocalele();
const int loc_dof_v = dof_mat_v * nLocBas;
const int loc_dof_p = dof_mat_p * nLocBas;

std::vector<std::vector<double>> array_cur_velo_sols(tm_RK_ptr->get_RK_step());
std::vector<std::vector<double>> array_cur_pres_sols(tm_RK_ptr->get_RK_step());

for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
{
  array_cur_velo_sols[ii] = cur_velo_sols[ii] -> GetLocalArray();
  array_cur_pres_sols[ii] = cur_pres_sols[ii] -> GetLocalArray();
}

const std::vector<double> array_cur_dot_velo = cur_dot_velo -> GetLocalArray();
const std::vector<double> array_cur_velo = cur_velo -> GetLocalArray();
const std::vector<double> array_pre_velo = pre_velo -> GetLocalArray();
const std::vector<double> array_cur_pres = cur_pres -> GetLocalArray();

double * ectrl_x = new double [nLocBas];
double * ectrl_y = new double [nLocBas];
double * ectrl_z = new double [nLocBas];
PetscInt * row_idx_v = new PetscInt [nLocBas * dof_mat_v];
PetscInt * row_idx_p = new PetscInt [nLocBas * dof_mat_p];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = locien->get_LIEN(ee);

  std::vector<std::vector<double>> local_cur_velo_sols(tm_RK_ptr->get_RK_step());    
  std::vector<std::vector<double>> local_cur_pres_sols(tm_RK_ptr->get_RK_step()); 

  for(int ii = 0; ii < tm_RK_ptr->get_RK_step(); ++ii)
  {
    local_cur_velo_sols[ii] = GetLocal( array_cur_velo_sols[ii], IEN_e, nLocBas, 3 );
    local_cur_pres_sols[ii] = GetLocal( array_cur_pres_sols[ii], IEN_e, nLocBas, 1 );
  }

  const std::vector<double> local_cur_dot_velo = GetLocal( array_cur_dot_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_cur_velo = GetLocal( array_cur_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_pre_velo = GetLocal( array_pre_velo, IEN_e, nLocBas, 3 );
  const std::vector<double> local_cur_pres = GetLocal( array_cur_pres, IEN_e, nLocBas, 1 );

  fnode->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  locassem->Assem_Tangent_Residual_Pressure(curr_time, dt, tm_RK_ptr, local_cur_dot_velo, 
      local_cur_velo_sols, local_cur_velo, local_cur_pres_sols, local_pre_velo, 
      local_cur_pres, ectrl_x, ectrl_y, ectrl_z);

      for(int ii=0; ii<nLocBas; ++ii)
      {
        const int loc_index  = locien->get_LIEN(ee, ii);
  
        row_idx_v[3*ii]   = 3 * nbc->get_LID( 1, loc_index );
        row_idx_v[3*ii+1] = 3 * nbc->get_LID( 2, loc_index ) + 1;
        row_idx_v[3*ii+2] = 3 * nbc->get_LID( 3, loc_index ) + 2;
  
        row_idx_p[ii] = nbc->get_LID( 0, loc_index );
      }

      MatSetValues(subK[0], loc_dof_p, row_idx_p, loc_dof_p, row_idx_p, locassem->Tangent0, ADD_VALUES);
      MatSetValues(subK[1], loc_dof_p, row_idx_p, loc_dof_v, row_idx_v, locassem->Tangent1, ADD_VALUES);
      MatSetValues(subK[2], loc_dof_v, row_idx_v, loc_dof_p, row_idx_p, locassem->Tangent2, ADD_VALUES);
      MatSetValues(subK[3], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent3, ADD_VALUES);
      MatSetValues(subK[4], loc_dof_v, row_idx_v, loc_dof_v, row_idx_v, locassem->Tangent4, ADD_VALUES);
      
      VecSetValues(subG[0], loc_dof_p, row_idx_p, locassem->Residual0, ADD_VALUES);
      VecSetValues(subG[1], loc_dof_v, row_idx_v, locassem->Residual1, ADD_VALUES);
}

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_idx_v; row_idx_v = nullptr;
delete [] row_idx_p; row_idx_p = nullptr;

NatBC_G_HERK_Pressure( curr_time, dt );

VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);

EssBC_KG();

MatAssemblyBegin(subK[0], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[0], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[1], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[1], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[2], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[2], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[3], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[3], MAT_FINAL_ASSEMBLY);
MatAssemblyBegin(subK[4], MAT_FINAL_ASSEMBLY); MatAssemblyEnd(subK[4], MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(subG[0]); VecAssemblyEnd(subG[0]);
VecAssemblyBegin(subG[1]); VecAssemblyEnd(subG[1]);
}

void PGAssem_Block_NS_FEM_HERK::NatBC_G_HERK_Sub( const double &curr_time, const double &dt,
  const int &substep_index,
  const ITimeMethod_RungeKutta * const &tm_RK_ptr )
{
int * LSIEN = new int [snLocBas];
double * sctrl_x = new double [snLocBas];
double * sctrl_y = new double [snLocBas];
double * sctrl_z = new double [snLocBas];
PetscInt * srow_index_v = new PetscInt [snLocBas * dof_mat_v];

for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
{
  const int num_sele = ebc -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc -> get_SIEN(ebc_id, ee, LSIEN);

    ebc -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    locassem->Assem_Residual_EBC_HERK_Sub(ebc_id, curr_time, dt, substep_index, tm_RK_ptr,
        sctrl_x, sctrl_y, sctrl_z);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat_v; ++mm)
        srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc -> get_LID(mm+1, LSIEN[ii]) + mm;
    }

    VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, locassem->Residual1, ADD_VALUES);
  }
}

delete [] LSIEN; LSIEN = nullptr;
delete [] sctrl_x; sctrl_x = nullptr;
delete [] sctrl_y; sctrl_y = nullptr;
delete [] sctrl_z; sctrl_z = nullptr;
delete [] srow_index_v; srow_index_v = nullptr;
}

void PGAssem_Block_NS_FEM_HERK::NatBC_G_HERK_Final( const double &curr_time, const double &dt,
  const ITimeMethod_RungeKutta * const &tm_RK_ptr )
{
int * LSIEN = new int [snLocBas];
double * sctrl_x = new double [snLocBas];
double * sctrl_y = new double [snLocBas];
double * sctrl_z = new double [snLocBas];
PetscInt * srow_index_v = new PetscInt [snLocBas * dof_mat_v];

for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
{
  const int num_sele = ebc -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc -> get_SIEN(ebc_id, ee, LSIEN);

    ebc -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    locassem->Assem_Residual_EBC_HERK_Final(ebc_id, curr_time, dt, tm_RK_ptr,
        sctrl_x, sctrl_y, sctrl_z);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat_v; ++mm)
        srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc -> get_LID(mm+1, LSIEN[ii]) + mm;
    }

    VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, locassem->sur_Residual1, ADD_VALUES);
  }
}

delete [] LSIEN; LSIEN = nullptr;
delete [] sctrl_x; sctrl_x = nullptr;
delete [] sctrl_y; sctrl_y = nullptr;
delete [] sctrl_z; sctrl_z = nullptr;
delete [] srow_index_v; srow_index_v = nullptr;
}

void PGAssem_Block_NS_FEM_HERK::NatBC_G_HERK_Pressure( const double &curr_time, const double &dt )
{
int * LSIEN = new int [snLocBas];
double * sctrl_x = new double [snLocBas];
double * sctrl_y = new double [snLocBas];
double * sctrl_z = new double [snLocBas];
PetscInt * srow_index_v = new PetscInt [snLocBas * dof_mat_v];

for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
{
  const int num_sele = ebc -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc -> get_SIEN(ebc_id, ee, LSIEN);

    ebc -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    locassem->Assem_Residual_EBC_HERK_Pressure(ebc_id, curr_time, dt,
        sctrl_x, sctrl_y, sctrl_z);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat_v; ++mm)
        srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc -> get_LID(mm+1, LSIEN[ii]) + mm;
    }

    VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, locassem->Residual1, ADD_VALUES);
  }
}

delete [] LSIEN; LSIEN = nullptr;
delete [] sctrl_x; sctrl_x = nullptr;
delete [] sctrl_y; sctrl_y = nullptr;
delete [] sctrl_z; sctrl_z = nullptr;
delete [] srow_index_v; srow_index_v = nullptr;
}

// EOF
