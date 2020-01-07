#include "PGAssem_FSI_FEM.hpp"

PGAssem_FSI_FEM::PGAssem_FSI_FEM( 
    IPLocAssem * const &locassem_f_ptr,
    IPLocAssem * const &locassem_s_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    IAGlobal_Mesh_Info const * const &agmi_ptr,
    ALocal_Elem const * const &alelem_ptr,
    ALocal_IEN const * const &aien_ptr,
    APart_Node const * const &pnode_ptr,
    ALocal_NodalBC const * const &part_nbc,
    ALocal_EBC const * const &part_ebc )
{
  nLocBas = agmi_ptr->get_nLocBas();
  dof_sol = pnode_ptr->get_dof(); // pnode_ptr stores dofNum
  dof_mat = locassem_f_ptr->get_dof_mat(); // locassem_ptr defines the matrix size
  num_ebc = part_ebc->get_num_ebc();

  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_sol != locassem_f_ptr->get_dof(),
      "PGAssem_FSI_FEM::dof_sol != locassem_f_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_sol != locassem_s_ptr->get_dof(),
      "PGAssem_FSI_FEM::dof_sol != locassem_f_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != part_nbc->get_dofMat(),
      "PGAssem_FSI_FEM::dof_mat != part_nbc->get_dofMat(). \n");

  SYS_T::print_fatal_if(dof_mat !=  locassem_s_ptr->get_dof_mat() ,
      "PGAssem_FSI_FEM::dof_mat != locassem_s_ptr->get_dof_mat. \n");

  // The users shall make sure that the surface element's number 
  // of local basis are the same. Here we take the value from the 
  // first ebc face.
  // This is an assumption in this assembly routine.
  if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  else snLocBas = 0;
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_FSI_FEM, snLocBas has to be uniform. \n");
  }

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlgn       = pnode_ptr->get_nlocghonode();
  const int nlocrow    = dof_mat * nlocalnode;

  // Allocate the AIJ matrix
  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];
  const int empirical_neibor_number = 250;
  Get_dnz_onz( nlocalnode,  empirical_neibor_number,
      part_nbc, dnnz, onnz );

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, dnnz, 0, onnz, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  delete [] dnnz; dnnz = NULL;
  delete [] onnz; onnz = NULL;

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");
  Release_nonzero_err_str();

  row_index = new PetscInt [nLocBas * dof_mat];

  array_a = new double [nlgn * dof_sol];
  array_b = new double [nlgn * dof_sol];

  local_a = new double [nLocBas * dof_sol];
  local_b = new double [nLocBas * dof_sol];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];

  if(num_ebc > 0)
  {
    LSIEN = new int [snLocBas];
    local_as = new double [dof_sol * snLocBas];
    local_bs = new double [dof_sol * snLocBas];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index = new PetscInt [dof_mat * snLocBas];
  }

  // Initialize array_b for nonzero assembly
  for(int ii=0; ii<nlgn*dof_sol; ++ii) array_b[ii] = 0.0;

  // Now we run a nonzero estimate trial assembly
  Assem_nonzero_estimate( alelem_ptr, locassem_f_ptr, locassem_s_ptr, 
      elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ebc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation
  
  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}


PGAssem_FSI_FEM::~PGAssem_FSI_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
  delete [] row_index; row_index = NULL;
  delete [] array_a;     array_a = NULL;
  delete [] array_b;     array_b = NULL;
  delete [] local_a;     local_a = NULL;
  delete [] local_b;     local_b = NULL;
  delete [] IEN_e;       IEN_e = NULL;
  delete [] ectrl_x;     ectrl_x = NULL;
  delete [] ectrl_y;     ectrl_y = NULL;
  delete [] ectrl_z;     ectrl_z = NULL;

  if(num_ebc > 0)
  {
    delete [] LSIEN; LSIEN = NULL;
    delete [] local_as; local_as = NULL;
    delete [] local_bs; local_bs = NULL;
    delete [] sctrl_x; sctrl_x = NULL;
    delete [] sctrl_y; sctrl_y = NULL;
    delete [] sctrl_z; sctrl_z = NULL;
    delete [] srow_index; srow_index = NULL;
  }
}


void PGAssem_FSI_FEM::Get_dnz_onz( const int &nlocnode,
    const int &empirical_neighbor_node_number,
    const ALocal_NodalBC * const &nbc_ptr,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nzbase = dof_mat * empirical_neighbor_node_number;

  Vec vdnz, vonz;
  VecCreateMPI(PETSC_COMM_WORLD, dof_mat * nlocnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof_mat * nlocnode, PETSC_DETERMINE, &vonz);
  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  VecSetOption(vdnz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(vonz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  int row;
  for(int ii=0; ii<nlocnode; ++ii)
  {
    for(int mm=0; mm<dof_mat; ++mm)
    {
      row = nbc_ptr->get_LID(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  for(int mm=0; mm<dof_mat; ++mm)
  {
    // Check the master nodes for each d.o.f.
    const int num_master = nbc_ptr->get_Num_LPM(mm);
    for(int ii=0; ii<num_master; ++ii)
    {
      row = nbc_ptr->get_LocalMaster(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);


  for(int mm=0; mm<dof_mat; ++mm)
  {
    // Check the Dirichlet nodes for each d.o.f.
    const int num_dir = nbc_ptr->get_Num_LD(mm);
    for(int ii=0; ii<num_dir; ++ii)
    {
      row = nbc_ptr->get_LDN(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
      VecSetValue(vonz, row, 0.0, INSERT_VALUES);
    }

    // Check the slave nodes for each d.o.f.
    const int num_slave = nbc_ptr->get_Num_LPS(mm);
    for(int ii=0; ii<num_slave; ++ii)
    {
      row = nbc_ptr->get_LPSN(mm, ii) * dof_mat + mm;
      VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
      VecSetValue(vonz, row, 2.0, INSERT_VALUES);
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  PetscInt mat_length;
  VecGetSize(vdnz, &mat_length);

  const int max_dnz = dof_mat * nlocnode;
  const int max_onz = mat_length - dof_mat * nlocnode;

  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof_mat*nlocnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    if(dnz[ii] > max_dnz)
      dnz[ii] = max_dnz;
  }
  VecRestoreArray(vdnz, &array_d);

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof_mat*nlocnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);
    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}

void PGAssem_FSI_FEM::EssBC_KG( 
    const ALocal_NodalBC * const &nbc_part, const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);

  if(local_dir > 0)
  {
    int row, col;
    for(int i=0; i<local_dir; ++i)
    {
      row = nbc_part->get_LDN(field, i) * dof_mat + field;
      col = row;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(field);
  if(local_sla > 0)
  {
    int row, col;
    for(int i=0; i<local_sla; ++i)
    {
      row = nbc_part->get_LPSN(field, i) * dof_mat + field;
      col = nbc_part->get_LPMN(field, i) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_FSI_FEM::EssBC_G( 
    const ALocal_NodalBC * const &nbc_part, const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);
  const int local_sla = nbc_part->get_Num_LPS(field);
  int row;
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      row = nbc_part->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      row = nbc_part->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_FSI_FEM::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_f_ptr,
    IPLocAssem * const &lassem_s_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  int loc_index, lrow_index, offset1;

  lassem_f_ptr->Assem_Estimate();
  lassem_s_ptr->Assem_Estimate();

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      loc_index  = lien_ptr->get_LIEN(e, i);

      offset1 = dof_mat * i;

      for(int m=0; m<dof_mat; ++m)
      {
        lrow_index = nbc_part->get_LID( m, loc_index );

        row_index[offset1 + m] = dof_mat * lrow_index + m;
      }
    }

    // If the element is fluid, call fluid local assembly;
    // otherwise, call solid local assembly
    if(alelem_ptr->get_elem_tag(e) == 0)
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_f_ptr->Tangent, ADD_VALUES);
    }
    else
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_s_ptr->Tangent, ADD_VALUES);
    }
  }

  // Resis BC for K and G
  // choose an arbitrary (0.1) time step size just to get the nonzero pattern
  NatBC_Resis_KG(0.1, lassem_f_ptr, elements, quad_s, nbc_part, ebc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  
  for(int fie=0; fie<dof_mat; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_FSI_FEM::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_f_ptr,
    IPLocAssem * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );

  // array_b here stores the sol_a's local copy. It will be used
  // in the NatBC_G calculation.
  sol_a->GetLocalArray( array_b, node_ptr );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // If elem_tag is 0, do fluid assembly, else do solid
    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      lassem_f_ptr->Assem_Mass_Residual( local_a, elementv,
          ectrl_x, ectrl_y, ectrl_z, quad_v );
    }
    else
    {
      lassem_s_ptr->Assem_Mass_Residual( local_a, elementv,
          ectrl_x, ectrl_y, ectrl_z, quad_v );
    }

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof_mat * ii;

      for(int mm=0; mm<dof_mat; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);
        row_index[offset1+mm] = dof_mat * lrow_index + mm;
      }
    }

    // If elem_tag is 0, use fluid pointer, else use solid pointer 
    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_f_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_f_ptr->Residual, ADD_VALUES);
    }
    else
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_s_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_s_ptr->Residual, ADD_VALUES);
    }
  }

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof_mat; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_FSI_FEM::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_f_ptr,
    IPLocAssem * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // If elem tag is zero, do fluid, otherwise, do solid
    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      lassem_f_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }
    else
    {
      lassem_s_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof_mat * ii;
      for(int mm=0; mm<dof_mat; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);
        row_index[offset1+mm] = dof_mat * lrow_index + mm;
      }
    }

    // If elem tag is zero, use fluid, otherwise, use solid
    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      VecSetValues(G, loc_dof, row_index, lassem_f_ptr->Residual, ADD_VALUES);
    }
    else
    {
      VecSetValues(G, loc_dof, row_index, lassem_s_ptr->Residual, ADD_VALUES);
    }
  }

  // Resistance BC for G
  NatBC_Resis_G( lassem_f_ptr, elements, quad_s, nbc_part, ebc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof_mat; ++fie) EssBC_G( nbc_part, fie );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_FSI_FEM::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_f_ptr,
    IPLocAssem * const &lassem_s_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    // If elem tag is zero, do fluid assembly; else do solid
    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      lassem_f_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }
    else
    {
      lassem_s_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
          elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);
    }

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof_mat * ii;

      for(int mm=0; mm<dof_mat; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);

        row_index[offset1 + mm] = dof_mat * lrow_index + mm;
      }
    }

    if(alelem_ptr->get_elem_tag(ee) == 0)
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_f_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_f_ptr->Residual, ADD_VALUES);
    }
    else
    {
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_s_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_s_ptr->Residual, ADD_VALUES);
    }
  }

  // Resistance BC for K and G
  NatBC_Resis_KG( dt, lassem_f_ptr, elements, quad_s, nbc_part, ebc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof_mat; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_FSI_FEM::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const int &in_loc_dof,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const lien_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part )
{
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      lassem_f_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          local_as, local_bs, element_s, sctrl_x, sctrl_y, sctrl_z,
          quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        int loc_index = LSIEN[ii];
        int offset1 = dof_mat * ii;
        for(int mm=0; mm<dof_mat; ++mm)
        {
          int lrow_index = nbc_part -> get_LID(mm, loc_index);
          srow_index[offset1 + mm] = dof_mat * lrow_index + mm;
        }
      }
      VecSetValues(G, in_loc_dof, srow_index, lassem_f_ptr->Residual, ADD_VALUES);
    }
  }
}


double PGAssem_FSI_FEM::Get_flow_rate( IPLocAssem * const &lassem_f_ptr,
    FEAElement * const &ele_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id ) const
{
  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double sum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    GetLocal(array_b, LSIEN, snLocBas, local_bs);

    sum += lassem_f_ptr -> get_flowrate( local_bs, ele_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  return sum;
}


void PGAssem_FSI_FEM::NatBC_Resis_G(
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate flow rate for face with ebc_id
    const double flrate = Assem_surface_flowrate( sol, lassem_f_ptr,
        element_s, quad_s, node_ptr, ebc_part, ebc_id );

    // Get the pressure value on the outlet surfaces
    const double P_n = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, flrate );

    // P_n+alpha_f
    const double val = P_n + lassem_ptr->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);
      
      lassem_f_ptr->Assem_Residual_EBC_Resistance(ebc_id, val,
          local_bs, element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = lassem_f_ptr->Residual[4*ii+1];
        Res[3*ii+1] = lassem_f_ptr->Residual[4*ii+2];
        Res[3*ii+2] = lassem_f_ptr->Residual[4*ii+3];

        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1, LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2, LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3, LSIEN[ii]) + 3;
      }

      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }
  }

  delete [] Res; Res = NULL; delete [] srow_idx; srow_idx = NULL;
}

void PGAssem_FSI_FEM::NatBC_Resis_KG(
    const double &dt,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_f_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const double dd_dv = dt * lassem_f_ptr->get_model_para_1() * lassem_f_ptr->get_model_para_2();

  // Allocate the vector to hold the residual on each surface element
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];
  PetscScalar * Tan;
  PetscInt * scol_idx;
  double out_nx, out_ny, out_nz;
  std::vector<double> intNB;
  std::vector<int> map_Bj;

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate flow rate for face with ebc_id
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, node_ptr, ebc_part, ebc_id );

    // Get the pressure value on the outlet surface
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, flrate );

    // P_n+alpha_f
    const double resis_val = P_n + lassem_ptr->get_model_para_1() * (P_np1 - P_n);

    // Get m := dP/dQ
    const double m_val = gbc -> get_m( ebc_id, flrate );

    const int num_face_nodes = ebc_part -> get_num_face_nodes(ebc_id);
    if(num_face_nodes > 0)
    {
      Tan = new PetscScalar [snLocBas * 3 * num_face_nodes * 3];
      scol_idx = new PetscInt [num_face_nodes * 3];
      ebc_part -> get_outvec( ebc_id, out_nx, out_ny, out_nz );
      ebc_part -> get_intNA( ebc_id, intNB );
      ebc_part -> get_LID( ebc_id, map_Bj );
    }
    else
    {
      Tan = NULL;
      scol_idx = NULL;
    }

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      // For here, we scale the int_NA nx/y/z by factor 1
      lassem_f_ptr->Assem_Residual_EBC_Resistance(ebc_id, 1.0,
          local_bs, element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      // Residual vector is scaled by the resistance value
      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = resis_val * lassem_f_ptr->Residual[4*ii+1];
        Res[3*ii+1] = resis_val * lassem_f_ptr->Residual[4*ii+2];
        Res[3*ii+2] = resis_val * lassem_f_ptr->Residual[4*ii+3];
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            Tan[temp_row + 3*B + 0] = dd_dv * m_val * lassem_f_ptr->Residual[4*A+ii+1] * intNB[B] * out_nx;
            Tan[temp_row + 3*B + 1] = dd_dv * m_val * lassem_f_ptr->Residual[4*A+ii+1] * intNB[B] * out_ny;
            Tan[temp_row + 3*B + 2] = dd_dv * m_val * lassem_f_ptr->Residual[4*A+ii+1] * intNB[B] * out_nz;
          }
        }
      }

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1,LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2,LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3,LSIEN[ii]) + 3;
      }

      for(int ii=0; ii<num_face_nodes; ++ii)
      {
        scol_idx[ii*3+0] = dof_mat * map_Bj[ii*3+0] + 1;
        scol_idx[ii*3+1] = dof_mat * map_Bj[ii*3+1] + 2;
        scol_idx[ii*3+2] = dof_mat * map_Bj[ii*3+2] + 3;
      }

      MatSetValues(K, snLocBas*3, srow_idx, num_face_nodes*3, scol_idx, Tan, ADD_VALUES);
      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }

    if( num_face_nodes > 0 ) 
    {
      delete [] Tan; 
      Tan = NULL;
      delete [] scol_idx;
      scol_idx = NULL;
    }
  }

  delete [] Res; Res = NULL; delete [] srow_idx; srow_idx = NULL;
}


double PGAssem_FSI_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const APart_Node * const &pnode_ptr,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  const int nlgn = pnode_ptr -> get_nlocghonode();
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array, pnode_ptr );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN );

    ebc_part -> get_ctrlPts_xyz( ebc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    GetLocal( array, LSIEN, snLocBas, local );

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array; array = NULL;
  delete [] local; local = NULL;

  double sum = 0.0;

  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}


double PGAssem_FSI_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const APart_Node * const &pnode_ptr,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  const int nlgn = pnode_ptr->get_nlocghonode();
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array, pnode_ptr );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = NULL;
  delete [] local; local = NULL;

  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}

// EOF
