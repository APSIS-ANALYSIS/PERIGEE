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
      dof_mat_v*in_nz_estimate, NULL, &K00);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_v, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_v*in_nz_estimate, NULL, 
      dof_mat_p*in_nz_estimate, NULL, &K01);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_v, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL, 
      dof_mat_v*in_nz_estimate, NULL, &K10);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow_p, nlocrow_p, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat_p*in_nz_estimate, NULL, 
      dof_mat_p*in_nz_estimate, NULL, &K11);

  // Allocate the sub-vectors
  VecCreate(PETSC_COMM_WORLD, &G0);
  VecCreate(PETSC_COMM_WORLD, &G1);

  VecSetSizes(G0, nlocrow_v, PETSC_DECIDE);
  VecSetSizes(G1, nlocrow_p, PETSC_DECIDE);

  VecSetFromOptions(G0);
  VecSetFromOptions(G1);

  VecSet(G0, 0.0);
  VecSet(G1, 0.0);
  
  VecSetOption(G0, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(G1, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);



}


PGAssem_2x2Block_NS_FEM::~PGAssem_2x2Block_NS_FEM()
{}


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
