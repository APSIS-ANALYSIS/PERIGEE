#include "PGAssem_Mesh.hpp"

PGAssem_Mesh::PGAssem_Mesh( IPLocAssem * const &locassem_ptr,
    IAGlobal_Mesh_Info const * const &agmi_ptr,
    ALocal_Elem const * const &alelem_ptr,
    ALocal_IEN const * const &aien_ptr,
    APart_Node const * const &pnode_ptr,
    ALocal_NodalBC const * const &part_nbc,
    ALocal_EBC const * const &part_ebc )
{}


PGAssem_Mesh::~PGAssem_Mesh()
{}


void PGAssem_Mesh::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part )
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

  for(int mm=0; mm<dof; ++mm) EssBC_KG( nbc_part, mm );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Mesh::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{}


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
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{}


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
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{}


void PGAssem_Mesh::EssBC_KG(
    const ALocal_NodalBC * const &nbc_part, 
    const int &field )
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


void PGAssem_Mesh::EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field )
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


void PGAssem_Mesh::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof * snLocBas];

  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
      
      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
        for(int mm=0; mm<dof; ++mm)
          srow_index[dof * ii + mm] =  nbc_part -> get_LID(mm, LSIEN[ii]);
  
      VecSetValues(G, dof*snLocBas, srow_index, lassem_ptr->Residual, ADD_VALUES);  
    }
  }

  delete [] LSIEN;      LSIEN = nullptr;
  delete [] sctrl_x;    sctrl_x = nullptr;
  delete [] sctrl_y;    sctrl_y = nullptr;
  delete [] sctrl_z;    sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}

// EOF
