#include "PGAssem_HERK_Block_NS_FEM.hpp"

PGAssem_HERK_Block_NS_FEM::PGAssem_HERK_Block_NS_FEM(
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    std::unique_ptr<ALocal_NBC> in_nbc,
    std::unique_ptr<ALocal_EBC> in_ebc,
    std::unique_ptr<ALocal_WeakBC> in_wbc,
    std::unique_ptr<IPLocAssem_2x2Block> in_locassem,    
    const int &in_nz_estimate )
: locien( std::move(in_locien) ),
  locelem( std::move(in_locelem) ),
  fnode( std::move(in_fnode) ),
  pnode( std::move(in_pnode) ),
  nbc( std::move(in_nbc) ),
  ebc( std::move(in_ebc) ),
  wbc( std::move(in_wbc) ),
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
      "PGAssem_NS_FEM::dof_sol != locassem->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat_v + dof_mat_p != nbc->get_dof_LID(),
      "PGAssem_NS_FEM::dof_mat != nbc->get_dof_LID(). \n");
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_NS_FEM, snLocBas has to be uniform. \n");
  }

  const int nlocrow_v  = dof_mat_v * pnode->get_nlocalnode(); 
  const int nlocrow_p  = dof_mat_p * pnode->get_nlocalnode(););

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

  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_HERK_Block_NS_FEM::~PGAssem_HERK_Block_NS_FEM()
{
  VecDestroy(&subG[0]); VecDestroy(&subG[1]); 
  MatDestroy(&subK[0]); MatDestroy(&subK[1]);
  MatDestroy(&subK[2]); MatDestroy(&subK[3]);
  MatDestroy(&subK[4]);
}

void PGAssem_HERK_Block_NS_FEM::EssBC_KG( const int &field )
{
  const int local_dir = nbc->get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc->get_LDN(field, i) * dof_mat + field;
      
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(field);
  if(local_sla > 0)
  {
    for(int i=0; i<local_sla; ++i)
    {
      const int row = nbc->get_LPSN(field, i) * dof_mat + field;
      const int col = nbc->get_LPMN(field, i) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_NS_FEM::EssBC_G( const int &field )
{
  const int local_dir = nbc->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc->get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_HERK_Block_NS_FEM::Assem_nonzero_estimate()
{
  const int nElem = locelem->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  locassem->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      const int loc_index  = locien->get_LIEN(e, i);

      for(int m=0; m<dof_mat; ++m)
        row_index[dof_mat * i + m] = dof_mat * nbc->get_LID( m, loc_index ) + m;
    }
    
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        locassem->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  // Create a temporary zero solution vector to feed Natbc_Resis_KG
  PDNSolution * temp = new PDNSolution_NS( pnode.get(), 0, false );

  // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  NatBC_Resis_KG( 0.0, 0.1, temp, temp, gbc );

  delete temp;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_HERK_Block_NS_FEM::Assem_mass_residual(
    const PDNSolution * const &sol_a)
{
  const int nElem = locelem->get_nlocalele();
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
    locien->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Mass_Residual( local_a, ectrl_x, ectrl_y, ectrl_z );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc -> get_LID(mm, IEN_e[ii]) + mm;
    }
    
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        locassem->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, locassem->Residual, ADD_VALUES);
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
  Weak_EssBC_G(0, 0, sol_a);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_HERK_Block_NS_FEM::Assem_tangent_residual_substep(
  const int &substep_index,
  PDNSolution ** const &cur_velo_sols,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution ** const &pre_velo_sols,
  PDNSolution * const &pre_velo,
  PDNSolution ** const &pre_pres_sols,
  PDNSolution * const &pre_velo_before,    
  const Runge_Kutta_Butcher * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = alelem_ptr->get_nlocalele();
const int loc_dof = dof_mat * nLocBas;

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
PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = lien_ptr->get_LIEN(ee);

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

  fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  lassem_ptr->Assem_Tangent_Residual_Substep(curr_time, dt, substep_index, tm_RK_ptr, local_cur_velo_sols, local_cur_pres_sols,
      local_pre_velo_sols, local_pre_pres_sols, local_pre_velo, local_pre_velo_before,
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

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_index; row_index = nullptr;

// Backflow stabilization residual & tangent contribution
// BackFlow_KG( dt, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

// Resistance type boundary condition
// NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
//    nbc_part, ebc_part, gbc );

NatBC_G_HERK_Sub( curr_time, dt, substep_index, tm_RK_ptr, lassem_ptr, elements, quad_s, 
   nbc_part, ebc_part );

// Weakly enforced no-slip boundary condition
// If wall_model_type = 0, it will do nothing.
// Weak_EssBC_KG(curr_time, dt, sol_b, lassem_ptr, elementvs, quad_s,
//   lien_ptr, fnode_ptr, nbc_part, wbc_part);

VecAssemblyBegin(G);
VecAssemblyEnd(G);

for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(G);
VecAssemblyEnd(G);
}

void PGAssem_HERK_Block_NS_FEM::Assem_tangent_residual_finalstep(
  PDNSolution ** const &cur_velo_sols,
  PDNSolution * const &cur_velo,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution ** const &pre_velo_sols,
  PDNSolution * const &pre_velo,
  PDNSolution ** const &pre_pres_sols,
  PDNSolution * const &pre_velo_before,    
  const Runge_Kutta_Butcher * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = alelem_ptr->get_nlocalele();
const int loc_dof = dof_mat * nLocBas;

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
PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = lien_ptr->get_LIEN(ee);

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

  fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  lassem_ptr->Assem_Tangent_Residual_Laststep(curr_time, dt, tm_RK_ptr, local_cur_velo_sols, local_cur_velo,
      local_cur_pres_sols, local_pre_velo_sols, local_pre_velo, local_pre_pres_sols, local_pre_velo_before,
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

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_index; row_index = nullptr;

// Backflow stabilization residual & tangent contribution
// BackFlow_KG( dt, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

// Resistance type boundary condition
// NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
//    nbc_part, ebc_part, gbc );

NatBC_G_HERK_Last( curr_time, dt, tm_RK_ptr, lassem_ptr, elements, quad_s, 
   nbc_part, ebc_part );

// Weakly enforced no-slip boundary condition
// If wall_model_type = 0, it will do nothing.
// Weak_EssBC_KG(curr_time, dt, sol_b, lassem_ptr, elementvs, quad_s,
//   lien_ptr, fnode_ptr, nbc_part, wbc_part);

VecAssemblyBegin(G);
VecAssemblyEnd(G);

for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(G);
VecAssemblyEnd(G);
}

void PGAssem_HERK_Block_NS_FEM::Assem_tangent_residual_presstage(
  PDNSolution * const &cur_dot_velo,
  PDNSolution ** const &cur_velo_sols,
  PDNSolution * const &cur_velo,
  PDNSolution ** const &cur_pres_sols,
  PDNSolution * const &pre_velo,
  PDNSolution * const &cur_pres,    
  const Runge_Kutta_Butcher * const &tm_RK_ptr,
  const double &curr_time,
  const double &dt )
{
const int nElem = alelem_ptr->get_nlocalele();
const int loc_dof = dof_mat * nLocBas;

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
PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

for(int ee=0; ee<nElem; ++ee)
{
  const std::vector<int> IEN_e = lien_ptr->get_LIEN(ee);

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

  fnode_ptr->get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

  lassem_ptr->Assem_Tangent_Residual_Finalstep(curr_time, dt, tm_RK_ptr, local_cur_dot_velo, 
      local_cur_velo_sols, local_cur_velo, local_cur_pres_sols, local_pre_velo, 
      local_cur_pres, elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    for(int mm=0; mm<dof_mat; ++mm)
      row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_e[ii])+mm;
  }

  MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
      lassem_ptr->Tangent, ADD_VALUES);

  VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
}

delete [] ectrl_x; ectrl_x = nullptr;
delete [] ectrl_y; ectrl_y = nullptr;
delete [] ectrl_z; ectrl_z = nullptr;
delete [] row_index; row_index = nullptr;

// Backflow stabilization residual & tangent contribution
// BackFlow_KG( dt, sol_b, lassem_ptr, elements, quad_s, nbc_part, ebc_part );

// Resistance type boundary condition
// NatBC_Resis_KG( curr_time, dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, 
//    nbc_part, ebc_part, gbc );
NatBC_G_HERK_Final( curr_time, dt, lassem_ptr, elements, quad_s, 
   nbc_part, ebc_part );

// Weakly enforced no-slip boundary condition
// If wall_model_type = 0, it will do nothing.
// Weak_EssBC_KG(curr_time, dt, sol_b, lassem_ptr, elementvs, quad_s,
//   lien_ptr, fnode_ptr, nbc_part, wbc_part);

VecAssemblyBegin(G);
VecAssemblyEnd(G);

for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(G);
VecAssemblyEnd(G);
}

// EOF
