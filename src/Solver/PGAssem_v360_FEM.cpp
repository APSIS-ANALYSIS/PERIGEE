#include "PGAssem_v360_FEM.hpp"

PGAssem_v360_FEM::PGAssem_v360_FEM(
    IPLocAssem * const &locassem_ptr,
    IAGlobal_Mesh_Info const * const &agmi_ptr,
    ALocal_Elem const * const &alelem_ptr,
    ALocal_IEN const * const &aien_ptr,
    APart_Node const * const &pnode_ptr,
    ALocal_NodalBC const * const &part_nbc,
    ALocal_EBC const * const &part_ebc )
{
  nLocBas = agmi_ptr -> get_nLocBas();
  dof = locassem_ptr -> get_dof();
  
  if( dof != pnode_ptr->get_dof() )
  {
    SYS_T::commPrint("Warning: The dof from Model does not match the dof from preprocessor. \n");
    PetscPrintf(PETSC_COMM_WORLD, "         dof from Model is %d \n", dof);
    PetscPrintf(PETSC_COMM_WORLD, "         dof from preprocessor is %d \n", pnode_ptr->get_dof());
  }

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlgn       = pnode_ptr->get_nlocghonode();
  const int nlocrow    = dof * nlocalnode;

  // Allocate the AIJ matrix
  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];
  const int empirical_neibor_number = 160;
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
 
  // Allocate useful arrays
  row_index = new PetscInt [nLocBas * dof];

  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  local_a = new double [ dof * nLocBas ];
  local_b = new double [ dof * nLocBas ];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
  
  // Allocate surface integral used arrays
  num_ebc = part_ebc -> get_num_ebc();
  
  SYS_T::print_fatal_if(num_ebc != locassem_ptr->get_num_ebc_fun(),
      "Error: The number of ebc does not match with the number of functions implemented in the local assembly routine. \n");


  if(num_ebc > 0)
    snLocBas = part_ebc -> get_cell_nLocBas(0);
  else
    snLocBas = 0;

  // Make sure that the surface are comprised of the same type element,
  // at least number of local basis are the same. This is an assumption
  // in this global assembly.
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_v360_FEM, snLocBas has to be uniform. \n");

  if(num_ebc > 0)
  {
    LSIEN = new int [snLocBas];
    local_as = new double [ dof * snLocBas ];
    local_bs = new double [ dof * snLocBas ];
    sctrl_x = new double [snLocBas];
    sctrl_y = new double [snLocBas];
    sctrl_z = new double [snLocBas];
    srow_index = new PetscInt [snLocBas * dof];
  }

  // Now we run a nonzero estimate trial assembly
  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, aien_ptr,
     pnode_ptr, part_nbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);
  
  MatDestroy(&K); // Destroy the K with rough preallocation
 
  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}


PGAssem_v360_FEM::~PGAssem_v360_FEM()
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


void PGAssem_v360_FEM::Get_dnz_onz( const int &nElem,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nlocalnode = node_ptr->get_nlocalnode();
  const int nnode = node_ptr->get_nlocghonode();

  std::vector<int> numLocNode;
  
  hid_t file_id = H5Fopen("NumLocalNode.h", H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5reader = new HDF5_Reader( file_id );

  hid_t drank;
  hsize_t * ddims;
  int * intarray;

  h5reader->read_intArray( "/", "nln", drank, ddims, intarray );

  if(drank > 1){
    std::stringstream ss; ss<<"Error: "<<"nln"<<" at "<<"/"<<" "
        <<"NumLocalNode"<<".h5 is not an one-dimensional array. \n";
    SYS_T::commPrint(ss.str().c_str());
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  const int length = ddims[0];

  VEC_T::fillArray(numLocNode, intarray, length);

  delete [] ddims; delete [] intarray;
  delete h5reader;
  H5Fclose( file_id );

  std::vector<unsigned int> nlist;
  nlist.clear();
  nlist.resize(numLocNode.size()+1);
  nlist[0] = 0;
  for(unsigned int ii=1; ii<=numLocNode.size(); ++ii)
    nlist[ii] = nlist[ii-1] + numLocNode[ii-1];

  // numlocNode is not needed anymore. nlist will provide the nodal 
  // partition info.
  VEC_T::clean(numLocNode);

  // This vector stores each row's diagonal col index and off-diagonal col index
  std::vector<int> interfun_d, interfun_o;

  // This MPI vector stores globally collected diagonal and off-diagonal col
  // number
  Vec vdnz, vonz;

  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vonz);

  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  int row, col, ien_index, part_id_row, part_id_col;

  std::vector<int> elem4node;

  // loop for each row 
  for( int ii=0; ii<nnode; ++ii )
  {
    elem4node.clear();
    for(int ee=0; ee<nElem; ++ee)
    {
      if( lien_ptr->isNode_in_Elem(ee, ii) )
        elem4node.push_back(ee);
    }

    for( int mm=0; mm<dof; ++mm )
    {
      row = nbc_part->get_LID( mm, ii );
      interfun_d.clear();
      interfun_o.clear();

      // only calculate for nondirichlet node
      if(row >= 0)
      {
        // based on nlist, find the row node's corresponding processor id.
        part_id_row = Get_part_id(nlist, row);
        for(unsigned int ei=0; ei<elem4node.size(); ++ei)
        {
          int ee = elem4node[ei];
          for(int kk=0; kk<nLocBas; ++kk)
          {
            ien_index = lien_ptr->get_LIEN(ee,kk);
            for( int nn=0; nn<dof; ++nn )
            {
              col = nbc_part->get_LID( nn, ien_index );

              if( col>=0 )
              {
                // based on nlist, find the processor that the col node 
                // belong to. If they belong to the same processor, add 
                // to diagonal, otherwise, add to off-diagonal.
                part_id_col = Get_part_id(nlist, col);

                if( part_id_row == part_id_col )
                  interfun_d.push_back(dof*col + nn);
                else
                  interfun_o.push_back(dof*col + nn);
              }
            }
          }
        }
        // now the interfun_d and interfun_o have cached all the col that is
        // nonzero for this row, we sort them and remove repeated col number.
        VEC_T::sort_unique_resize(interfun_d);
        VEC_T::sort_unique_resize(interfun_o);

        // Finish calculating for each row, the real row index is
        // dof * row + mm
        VecSetValue(vdnz, dof*row+mm, double(interfun_d.size()), ADD_VALUES);
        VecSetValue(vonz, dof*row+mm, double(interfun_o.size()), ADD_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // We need to handle the Dirichlet and Periodic Slave nodes
  for(int mm=0; mm<dof; ++mm)
  {
    int local_dir = nbc_part->get_Num_LD(mm);
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = nbc_part->get_LDN(mm, ii) * dof + mm;
        VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
        VecSetValue(vonz, row, 0.0, INSERT_VALUES);
      }
    }

    int local_sla = nbc_part->get_Num_LPS(mm);

    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        int row = nbc_part->get_LPSN(mm,ii) * dof + mm;
        VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
        VecSetValue(vonz, row, 2.0, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // Get the global vector size
  PetscInt vec_size;
  VecGetSize(vdnz, &vec_size);

  // Now we get the globally collected dnz and onz number
  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    // the estimator above may overestimate for periodic master nodes.
    // if the number of nonzeros is greater than the dof * nlocalnode
    // reduce it to full diagonal rows. Otherwise PETSc will throw an
    // error message.
    if(dnz[ii] > dof*nlocalnode)
      dnz[ii] = dof * nlocalnode;                                   
  }
  VecRestoreArray(vdnz, &array_d);

  const int max_onz = vec_size - dof * nlocalnode;

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);

    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}


void PGAssem_v360_FEM::Get_dnz_onz( const int &nlocnode,
    const int &empirical_neighbor_node_number,
    const ALocal_NodalBC * const &nbc_ptr,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nzbase = dof * empirical_neighbor_node_number;

  Vec vdnz, vonz;
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocnode, PETSC_DETERMINE, &vonz);
  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  VecSetOption(vdnz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(vonz, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  int row;
  for(int ii=0; ii<nlocnode; ++ii)
  {
    for(int mm=0; mm<dof; ++mm)
    {
      row = nbc_ptr->get_LID(mm, ii) * dof + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  for(int mm=0; mm<dof; ++mm)
  {
    // Check the master nodes for each d.o.f.
    const int num_master = nbc_ptr->get_Num_LPM(mm);
    for(int ii=0; ii<num_master; ++ii)
    {
      row = nbc_ptr->get_LocalMaster(mm, ii) * dof + mm;
      VecSetValue(vdnz, row, double(nzbase), ADD_VALUES);
      VecSetValue(vonz, row, double(nzbase), ADD_VALUES);
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);


  for(int mm=0; mm<dof; ++mm)
  {
    // Check the Dirichlet nodes for each d.o.f.
    const int num_dir = nbc_ptr->get_Num_LD(mm);
    for(int ii=0; ii<num_dir; ++ii)
    {
      row = nbc_ptr->get_LDN(mm, ii) * dof + mm;
      VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
      VecSetValue(vonz, row, 0.0, INSERT_VALUES);
    }

    // Check the slave nodes for each d.o.f.
    const int num_slave = nbc_ptr->get_Num_LPS(mm);
    for(int ii=0; ii<num_slave; ++ii)
    {
      row = nbc_ptr->get_LPSN(mm, ii) * dof + mm;
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

  const int max_dnz = dof * nlocnode;
  const int max_onz = mat_length - dof * nlocnode;

  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof*nlocnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    if(dnz[ii] > max_dnz)
      dnz[ii] = max_dnz;
  }
  VecRestoreArray(vdnz, &array_d);

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof*nlocnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);
    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}


void PGAssem_v360_FEM::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  lassem_ptr->Assem_Estimate();

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      loc_index  = lien_ptr->get_LIEN(e, i);

      offset1 = dof * i;

      for(int m=0; m<dof; ++m)
      {
        lrow_index = nbc_part->get_LID( m, loc_index );

        row_index[offset1 + m] = dof * lrow_index + m;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  for( int fie=0; fie<dof; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_v360_FEM::EssBC_KG( 
    const ALocal_NodalBC * const &nbc_part, const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);

  // dirichlet condition
  if(local_dir > 0)
  {
    int row, col;
    for(int i=0; i<local_dir; ++i)
    {
      row = nbc_part->get_LDN(field, i) * dof + field;
      col = row;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
    }
  }
  
  // master-slave relation
  const int local_sla = nbc_part->get_Num_LPS(field);
  if(local_sla > 0)
  {
    int row, col;
    for(int i=0; i<local_sla; ++i)
    {
      row = nbc_part->get_LPSN(field, i) * dof + field;
      col = nbc_part->get_LPMN(field, i) * dof + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_v360_FEM::EssBC_G( 
    const ALocal_NodalBC * const &nbc_part, const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);
  const int local_sla = nbc_part->get_Num_LPS(field);
  int row;
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      row = nbc_part->get_LDN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      row = nbc_part->get_LPSN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_v360_FEM::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
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

      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          local_as, local_bs, element_s, sctrl_x, sctrl_y, sctrl_z,
          quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        int loc_index = LSIEN[ii];
        int offset1 = dof * ii;
        for(int mm=0; mm<dof; ++mm)
        {
          int lrow_index = nbc_part -> get_LID(mm, loc_index);
          srow_index[offset1 + mm] = dof * lrow_index + mm;
        }
      }

      VecSetValues(G, in_loc_dof, srow_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
}


void PGAssem_v360_FEM::Assem_mass_residual(
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
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  
  // array_b is used in NatBC_G call
  sol_a->GetLocalArray( array_b, node_ptr );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Mass_Residual( local_a, elementv,
        ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof * ii;

      for(int mm=0; mm<dof; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);
        row_index[offset1+mm] = dof * lrow_index + mm;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES); 
  }

  // Evaluate natural BC at time 0.0
  NatBC_G(0.0, 0.0, lassem_ptr, elements, dof*snLocBas, quad_s, lien_ptr,
      nbc_part, ebc_part);

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_v360_FEM::Assem_residual(
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
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof * ii;
      for(int mm=0; mm<dof; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);
        row_index[offset1+mm] = dof * lrow_index + mm;
      }
    }
    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  NatBC_G(curr_time, dt, lassem_ptr, elements, dof*snLocBas, quad_s,
      lien_ptr, nbc_part, ebc_part );
  
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie) EssBC_G( nbc_part, fie );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_v360_FEM::Assem_tangent_residual(
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
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN_e(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      loc_index = IEN_e[ii];
      offset1 = dof * ii;

      for(int mm=0; mm<dof; ++mm)
      {
        lrow_index = nbc_part -> get_LID(mm, loc_index);

        row_index[offset1 + mm] = dof * lrow_index + mm;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  NatBC_G(curr_time, dt, lassem_ptr, elements, dof*snLocBas, quad_s,
      lien_ptr, nbc_part, ebc_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie) EssBC_KG( nbc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


//EOF
