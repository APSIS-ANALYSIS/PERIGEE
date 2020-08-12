#include "PGAssem_2x2Block_NS_FEM.hpp"

PGAssem_2x2Block_NS_FEM::PGAssem_2x2Block_NS_FEM(
    IPLocAssem_2x2Block * const &locassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &part_nbc,
    const ALocal_EBC * const &part_ebc,
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat_v( locassem_ptr->get_dof_mat_0() ), 
  dof_mat_p( locassem_ptr->get_dof_mat_1() ),
  num_ebc( part_ebc->get_num_ebc() )
{
  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),
      "PGAssem_NS_FEM::dof_sol != locassem_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat_v + dof_mat_p != part_nbc->get_dofMat(),
      "PGAssem_NS_FEM::dof_mat != part_nbc->get_dofMat(). \n");

  // Make sure that the surface element's number of local basis
  // are the same by the users.
  if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  else snLocBas = 0;

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_NS_FEM, snLocBas has to be uniform. \n");
  }

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlgn       = pnode_ptr->get_nlocghonode();
  const int nlocrow_v  = dof_mat_v * nlocalnode;
  const int nlocrow_p  = dof_mat_p * nlocalnode;

  // Allocate the block matrices
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
      dof_mat_v*in_nz_estimate, NULL, &subK[0]);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
      dof_mat_p*in_nz_estimate, NULL, &subK[1]);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL, 
      dof_mat_v*in_nz_estimate, NULL, &subK[2]);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL, 
      dof_mat_p*in_nz_estimate, NULL, &subK[3]);

  // Allocate the sub-vectors
  VecCreate(PETSC_COMM_WORLD, &subG[0]);
  VecCreate(PETSC_COMM_WORLD, &subG[1]);

  VecSetSizes(subG[0], nlocrow_v, PETSC_DECIDE);
  VecSetSizes(subG[1], nlocrow_p, PETSC_DECIDE);

  VecSetFromOptions(subG[0]);
  VecSetFromOptions(subG[1]);

  VecSet(subG[0], 0.0);
  VecSet(subG[1], 0.0);
  
  VecSetOption(subG[0], VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(subG[1], VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  // Temporarily ignore new entry allocation
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  // Allocate containers
  row_index_v = new PetscInt [nLocBas * dof_mat_v];
  row_index_p = new PetscInt [nLocBas * dof_mat_p];

  array_v = new double [nlgn * dof_mat_v];
  array_p = new double [nlgn * dof_mat_p];
  array_dot_v = new double [nlgn * dof_mat_v];
  array_dot_p = new double [nlgn * dof_mat_p];

  local_v = new double [nLocBas * dof_mat_v];
  local_p = new double [nLocBas * dof_mat_p];
  local_dot_v = new double [nLocBas * dof_mat_v];
  local_dot_p = new double [nLocBas * dof_mat_p];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];

  if(num_ebc > 0)
  {
    LSIEN = new int [snLocBas];
    local_vs = new double [snLocBas * dot_mat_v];
    local_ps = new double [snLocBas * dot_mat_p];
    local_dot_vs = new double [snLocBas * dot_mat_v];
    local_dot_ps = new double [snLocBas * dot_mat_p];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index_v = new PetscInt [snLocBas * dof_mat_v];
    srow_index_p = new PetscInt [snLocBas * dof_mat_p];
  }

  // Initialize values for array
  for(int ii=0; ii<nlgn*dof_mat_v; ++ii)
  {
    array_v[ii] = 0.0;
    array_dot_v[ii] = 0.0;
  }

  for(int ii=0; ii<nlgn*dof_mat_p; ++ii)
  {
    array_p[ii] = 0.0;
    array_dot_p[ii] = 0.0;
  }

  // Nonzero pattern assembly
  //Assem_nonzero_estimate( alelem_ptr, locassem_ptr,
  //    elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ebc, gbc );

}


PGAssem_2x2Block_NS_FEM::~PGAssem_2x2Block_NS_FEM()
{}


void PGAssem_2x2Block_NS_FEM::EssBC_KG_v( const ALocal_NodalBC * const &nbc_part )
{
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc_part->get_Num_LD(field);

    if(local_dir > 0)
    {
      for(int i=0; i<local_dir; ++i)
      {
        const int row = nbc_part->get_LDN(field, i) * dof_mat_v + field - 1;
        VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
        MatSetValue(subK[0], row, row, 1.0, ADD_VALUES);
      }
    }

    const int local_sla = nbc_part->get_Num_LPS(field);
    if(local_sla > 0)
    {
      for(int i=0; i<local_sla; ++i)
      {
        const int row = nbc_part->get_LPSN(field, i) * dof_mat_v + field - 1;
        const int col = nbc_part->get_LPMN(field, i) * dof_mat_v + field - 1;
        MatSetValue(subK[0], row, col, 1.0, ADD_VALUES);
        MatSetValue(subK[0], row, row, -1.0, ADD_VALUES);
        VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_NS_FEM::EssBC_G_v( const ALocal_NodalBC * const &nbc_part )
{
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc_part->get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part->get_LDN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc_part->get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_part->get_LPSN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_NS_FEM::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{}










// EOF
