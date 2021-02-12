#include "PGAssem_v360_noCache_NURBS.hpp"

PGAssem_v360_noCache_NURBS::PGAssem_v360_noCache_NURBS(
    const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const IALocal_BC * const &part_bc )
: nLocBas(agmi_ptr->get_nLocBas()) , dof(locassem_ptr->get_dof())
{
  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlocrow = dof * nlocalnode;
  const int nElem = alelem_ptr->get_nlocalele();

  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];

  SYS_T::commPrint("===> Estimate sparse nonzero structure. \n");
  Get_dnz_onz(nElem, aien_ptr, pnode_ptr, part_bc, dnnz, onnz);

  Init_petsc_360(dnnz, onnz, nlocrow);

  delete [] dnnz; delete [] onnz;

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");
  Release_nonzero_err_str();

  // allocate the frequently used arrays in global assembly
  const int nlgn = pnode_ptr->get_nlocghonode();
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  row_index = new PetscInt [dof * nLocBas];

  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
  ectrl_w = new double [nLocBas];
}



PGAssem_v360_noCache_NURBS::~PGAssem_v360_noCache_NURBS()
{
  VecDestroy(&G);
  MatDestroy(&K);
  delete [] row_index;
  delete [] array_a;
  delete [] array_b;
  delete [] local_a;
  delete [] local_b;
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] ectrl_w;
}


void PGAssem_v360_noCache_NURBS::Init_petsc_360(const PetscInt * const &dnz,
    const PetscInt * const &onz, const int &num_loc_row )
{
  SYS_T::commPrint("===> PETSc-3.6.0: MatCreateAIJ called. \n");
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, 0, dnz, 0, onz, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}



void PGAssem_v360_noCache_NURBS::Get_dnz_onz( const int &nElem,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const IALocal_BC * const &bc_part,
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

  if(drank > 1)
  {
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
      row = bc_part->get_LID( mm, ii );
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
              col = bc_part->get_LID( nn, ien_index );

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
    int local_dir = bc_part->get_Num_LD(mm);
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = bc_part->get_LDN(mm, ii) * dof + mm;
        VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
        VecSetValue(vonz, row, 0.0, INSERT_VALUES);
      }
    }

    int local_sla = bc_part->get_Num_LP(mm);

    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        int row = bc_part->get_LPSN(mm,ii) * dof + mm;
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



void PGAssem_v360_noCache_NURBS::EssBC_KG( const IALocal_BC * const &bc_part, 
    const int &field )
{
  // Check if dirichlet nodes exists within this partition
  // NOTE: we use ADD_VALUES here for matrix assembly, where we assumes that
  // the matrix is assemblyed with LID, which does nothing to the essential
  // boundary nodes, i.e., the boundary nodes' rows are zero rows.
  const int local_dir = bc_part->get_Num_LD(field);
  
  if(local_dir > 0)
  {
    int row, col;
    for(int i=0; i<local_dir; ++i)
    {
      row = bc_part->get_LDN(field, i) * dof + field;
      col = row;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
    }
  }

  // Check if periodic slave nodes exists in this partition
  const int local_sla = bc_part->get_Num_LP(field);
  if(local_sla > 0)
  {
    int row, col;
    for(int i=0; i<local_sla; ++i)
    {
      row = bc_part->get_LPSN(field, i) * dof + field;
      col = bc_part->get_LPMN(field, i) * dof + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}



void PGAssem_v360_noCache_NURBS::EssBC_G( const IALocal_BC * const &bc_part, 
    const int &field )
{
  const int local_dir = bc_part->get_Num_LD(field);
  const int local_sla = bc_part->get_Num_LP(field);
  int row;
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      row = bc_part->get_LDN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      row = bc_part->get_LPSN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}



void PGAssem_v360_noCache_NURBS::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const IALocal_BC * const &bc_part )
{  
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  // Generate a local matrix with 1.0 in all slots 
  lassem_ptr->Assem_Estimate();

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      loc_index  = lien_ptr->get_LIEN(e, i);
      
      offset1 = dof * i;

      for(int m=0; m<dof; ++m)
      {
        lrow_index = bc_part->get_LID( m, loc_index );

        row_index[offset1 + m] = dof * lrow_index + m;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  for( int fie=0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}



void PGAssem_v360_noCache_NURBS::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_BC * const &bc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  double ehx, ehy, ehz; // element mesh size in each direction

  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction

  sol_a->GetLocalArray( array_a, node_ptr );

  for(int ee=0; ee<nElem; ++ee)
  {
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a);

      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      lassem_ptr->Assem_Mass_Residual(local_a, ee, ehx, ehy, ehz, bs, bt, bu, 
          ectrl_x, ectrl_y, ectrl_z, ectrl_w, &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr); 

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];

        offset1 = dof * i; 

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);

          row_index[offset1 + m] = dof * lrow_index + m;
        }
      }
      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}



void PGAssem_v360_noCache_NURBS::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_BC * const &bc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  double ehx, ehy, ehz; // element mesh size in each direction

  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction

  for( int ee=0; ee<nElem; ++ee )
  {
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      // Obtain the mesh size in element ee
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      // Obtain the extraction operator in element ee
      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      // Obtain the geometrical infomation
      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      // Call local assembly routine to generate element matrix and vector
      lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, bs, bt, bu, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr);

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
  
        offset1 = dof * i;

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[offset1 + m] = dof * lrow_index + m;
        }
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
 
  // boundary integral on top 
  for(int ei=0; ei<bc_part->get_NumLE_Top(0); ++ei)
  {
    // get the element index whose top surface needs to perform integral
    int ee = bc_part->get_LTop_Elem(0, ei);
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      lassem_ptr->Assem_Residual_TopFace( curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0] );

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        offset1 = dof * i;
        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[offset1 + m] = dof * lrow_index + m;
        } 
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
  
  // boundary integral on bottom
  for(int ei=0; ei<bc_part->get_NumLE_Bot(0); ++ei)
  {
    // get the element index whose top surface needs to perform integral
    int ee = bc_part->get_LBottom_Elem(0, ei);
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      lassem_ptr->Assem_Residual_BotFace( curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0] );

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        offset1 = dof * i;
        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[offset1 + m] = dof * lrow_index + m;
        } 
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }


  // The Rest boundary integrals need to be completed

  // Finish the boundary integral assembly

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}





void PGAssem_v360_noCache_NURBS::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_BC * const &bc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index, offset1; // lcol_index;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );

  double ehx, ehy, ehz; // element mesh size in each direction

  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction


  // Volumetric integral
  for( int ee=0; ee<nElem; ++ee )
  {
    //rmeshsize = mSize->get_meshsize(ee);

    if( mSize->get_meshsize(ee) > 0.0 )
    {
      // Obtian the mesh size in element ee
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      // Obtain the extractor in element ee
      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      // Obtain the control points gemoetrical information
      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      // Call local assembly routine to generate element matrix and vector
      lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, bs, bt, bu, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr);


      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
         
        offset1 = dof * i;

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);

          row_index[offset1 + m] = dof * lrow_index + m;
        }
      }

      MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
          lassem_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }

  // Loop over top boundary face integral
  // Due to the design of the boundary element info, one has to specify the
  // location of boundary integral for the dof with index 0, even if the
  // boundary condition is applied for other dofs. In the global
  // assembly routine, I will only check the 0-index in the local bc object,
  // if it indicates that the bc element number is greater than 0, boundary
  // integral will be performed. This might be a bad design. However, if we
  // adopt weak-imposition of dirichlet bc in the future, this does not matter
  // at all. We will eventaully perform boundary integral for all dof's on all
  // faces for all boundary conditions. The difference of difference bc will
  // appear only in the local assembly routine.
  for(int ei=0; ei<bc_part->get_NumLE_Top(0); ++ei)
  {
    // get the element index whose top surface needs to perform integral
    int ee = bc_part->get_LTop_Elem(0, ei);
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      lassem_ptr->Assem_Residual_TopFace( curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0] );

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        
        offset1 = dof * i;

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[offset1 + m] = dof * lrow_index + m;
        } 
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
  
  for(int ei=0; ei<bc_part->get_NumLE_Bot(0); ++ei)
  {
    // get the element index whose top surface needs to perform integral
    int ee = bc_part->get_LBottom_Elem(0, ei);
    if( mSize->get_meshsize(ee) > 0.0 )
    {
      ehx = mSize->get_hx(ee);
      ehy = mSize->get_hy(ee);
      ehz = mSize->get_hz(ee);

      extractor->get_EXT_x(ee, ext_x);
      extractor->get_EXT_y(ee, ext_y);
      extractor->get_EXT_z(ee, ext_z);

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);

      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);

      lassem_ptr->Assem_Residual_BotFace( curr_time, dt, local_a, local_b,
          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
          &ext_x[0], &ext_y[0], &ext_z[0] );

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        offset1 = dof * i;
        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[offset1 + m] = dof * lrow_index + m;
        } 
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }


  // The Rest boundary integrals need to be completed

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  // Imposing essential boundary conditions
  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}




// EOF
