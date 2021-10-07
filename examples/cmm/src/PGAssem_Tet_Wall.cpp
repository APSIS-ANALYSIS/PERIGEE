#include "PGAssem_Tet_Wall.hpp"

PGAssem_Tet_Wall::PGAssem_Tet_Wall( 
    const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &part_nbc,
    const ALocal_EBC * const &part_ebc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat( locassem_ptr->get_dof_mat() ),
  nlgn( pnode_ptr->get_nlocghonode() ),
  snLocBas( part_ebc -> get_cell_nLocBas(0) )
{
  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),
      "PGAssem_Tet_CMM_GenAlpha::dof_sol != locassem_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != part_nbc->get_dofMat(),
      "PGAssem_Tet_CMM_GenAlpha::dof_mat != part_nbc->get_dofMat(). \n");

  const int nlocrow = dof_mat * pnode_ptr->get_nlocalnode();

  // Allocate the sparse matrix K
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat*in_nz_estimate, NULL, dof_mat*in_nz_estimate, NULL, &K);

  // Allocate the vector G
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate( alelem_ptr, aien_ptr, part_nbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation

  // Create Mat with precise preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}


PGAssem_Tet_Wall::~PGAssem_Tet_Wall()
{
  VecDestroy(&G);
  MatDestroy(&K);
}


void PGAssem_Tet_Wall::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_NodalBC * const &nbc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  PetscScalar * ones = new PetscScalar[4 * nLocBas * 4 * nLocBas];

  for(int ii=0; ii<16*nLocBas * nLocBas; ++ii) ones[ii] = 1.0;

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];
  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      const int loc_index  = lien_ptr->get_LIEN(e, i);

      for(int m=0; m<dof_mat; ++m)
        row_index[dof_mat * i + m] = dof_mat * nbc_part->get_LID( m, loc_index ) + m;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index, ones, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;
  delete [] ones;      ones      = nullptr;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_Wall::EssBC_KG( const ALocal_NodalBC * const &nbc_part, 
    const int &field )
{
  const int num_nbc = nbc_part->get_num_nbc();

  for(int nbc_id; nbc_id<num_nbc; ++nbc_id)
  {
    const int local_dir = nbc_part->get_Num_LD(nbc_id, field);

    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part->get_LDN(nbc_id, field, ii) * dof_mat + field;

        VecSetValue(G, row, 0.0, INSERT_VALUES);
        MatSetValue(K, row, row, 1.0, ADD_VALUES);
      }
    }
  }
}


void PGAssem_Tet_Wall::RingBC_KG(
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const int &dof, const int &nrow, const int &ncol,
    const PetscInt * const &row_index,
    const PetscInt * const &col_index,
    PetscScalar * const &Ke, 
    PetscScalar * const &Ge )
{
  const int ringbc_type = ringnbc_part -> get_ringbc_type();

  // Clamped rings
  if( ringbc_type == 0 ) {} 

  // Skew boundary conditions for in-plane motion of ring nodes
  else if( ringbc_type == 1 )
  {
    int pos = -1;  

    // Note: element tangent from NatBC_Resis_KG isn't a square matrix,
    //       ncol >= nrow.
    for( int ii = dof-1; ii < ncol; ii += dof )
    {    
      // Use velo-Z dof to determine ring nodes
      const int dnode = ( col_index[ii] - 3 ) / dof_mat;
      if( ringnbc_part->is_inLDN( dnode, pos ) )
      {    
        Matrix_3x3 Q = ringnbc_part->get_rotation_matrix( pos );
        Q.transpose(); // Skew-to-global transformation matrix

        for( int jj = dof-1; jj < nrow; jj += dof )
        {
          if( dnode != ( row_index[jj] - 3 ) / dof_mat )
          {
            Matrix_3x3 Ke_AB = Matrix_3x3(
              Ke[(jj-2)*ncol + (ii-2)], Ke[(jj-2)*ncol + (ii-1)], Ke[(jj-2)*ncol + ii], 
              Ke[(jj-1)*ncol + (ii-2)], Ke[(jj-1)*ncol + (ii-1)], Ke[(jj-1)*ncol + ii], 
              Ke[(jj-0)*ncol + (ii-2)], Ke[(jj-0)*ncol + (ii-1)], Ke[(jj-0)*ncol + ii]  );

            Ke_AB.MatMult(Ke_AB, Q);  // Ke_AB * Q

            // Update Ke
            Ke[(jj-2)*ncol + (ii-2)] = Ke_AB.xx(); Ke[(jj-2)*ncol + (ii-1)] = Ke_AB.xy(); Ke[(jj-2)*ncol + ii] = Ke_AB.xz(); 
            Ke[(jj-1)*ncol + (ii-2)] = Ke_AB.yx(); Ke[(jj-1)*ncol + (ii-1)] = Ke_AB.yy(); Ke[(jj-1)*ncol + ii] = Ke_AB.yz(); 
            Ke[(jj-0)*ncol + (ii-2)] = Ke_AB.zx(); Ke[(jj-0)*ncol + (ii-1)] = Ke_AB.zy(); Ke[(jj-0)*ncol + ii] = Ke_AB.zz(); 

            // Continuity eqn
            if( dof == dof_mat )
            {
              Vector_3 Ke_c = Vector_3( Ke[(jj-3)*ncol + (ii-2)], Ke[(jj-3)*ncol + (ii-1)], Ke[(jj-3)*ncol + ii] );
              Vector_3 rot_Ke_c = Q.VecMultT( Ke_c );  // rot_Ke_c = Ke_c^T * Q

              Ke[(jj-3)*ncol + (ii-2)] = rot_Ke_c.x(); Ke[(jj-3)*ncol + (ii-1)] = rot_Ke_c.y(); Ke[(jj-3)*ncol + ii] = rot_Ke_c.z(); 
            }
          }
        }
      }    
    }    

    for( int ii = dof-1; ii < nrow; ii += dof )
    {
      // Use velo-Z dof to determine ring nodes
      const int dnode = ( row_index[ii] - 3 ) / dof_mat;
      if( ringnbc_part->is_inLDN( dnode, pos ) )
      {
        // Global-to-skew transformation matrix
        Matrix_3x3 QT = ringnbc_part->get_rotation_matrix( pos );

        // Skew-to-global transformation matrix
        Matrix_3x3 Q( QT ); Q.transpose();

        for( int jj = dof-1; jj < ncol; jj += dof )
        {
          Matrix_3x3 Ke_AB = Matrix_3x3(
            Ke[(ii-2)*ncol + (jj-2)], Ke[(ii-2)*ncol + (jj-1)], Ke[(ii-2)*ncol + jj],
            Ke[(ii-1)*ncol + (jj-2)], Ke[(ii-1)*ncol + (jj-1)], Ke[(ii-1)*ncol + jj],
            Ke[(ii-0)*ncol + (jj-2)], Ke[(ii-0)*ncol + (jj-1)], Ke[(ii-0)*ncol + jj]  );

          if( dnode != ( col_index[jj] - 3 ) / dof_mat ) Ke_AB.MatMult( QT, Ke_AB );  // QT * Ke_AB
          else
          {
            Ke_AB.MatRot( Q );  // QT * Ke_AB * Q

            // Continuity eqn
            if( dof == dof_mat )
            {
              Vector_3 Ke_c = Vector_3( Ke[(ii-3)*ncol + (jj-2)], Ke[(ii-3)*ncol + (jj-1)], Ke[(ii-3)*ncol + jj] );
              Vector_3 rot_Ke_c = Q.VecMultT( Ke_c );  // rot_Ke_c = Ke_c^T * Q

              Ke[(ii-3)*ncol + (jj-2)] = rot_Ke_c.x(); Ke[(ii-3)*ncol + (jj-1)] = rot_Ke_c.y(); Ke[(ii-3)*ncol + jj] = rot_Ke_c.z();
            }
          }

          // Momentum eqns corresponding to pressure dof
          if( dof == dof_mat )
          {
            Vector_3 Ke_c = Vector_3( Ke[(ii-2)*ncol + (jj-3)], Ke[(ii-1)*ncol + (jj-3)], Ke[ii*ncol + (jj-3)] );
            Vector_3 rot_Ke_c = QT.VecMult( Ke_c );  // rot_Ke_c = Q^T * Ke_c

            Ke[(ii-2)*ncol + (jj-3)] = rot_Ke_c.x(); Ke[(ii-1)*ncol + (jj-3)] = rot_Ke_c.y(); Ke[ii*ncol + (jj-3)] = rot_Ke_c.z();
          }

          // Update Ke
          Ke[(ii-2)*ncol + (jj-2)] = Ke_AB.xx(); Ke[(ii-2)*ncol + (jj-1)] = Ke_AB.xy(); Ke[(ii-2)*ncol + jj] = Ke_AB.xz();
          Ke[(ii-1)*ncol + (jj-2)] = Ke_AB.yx(); Ke[(ii-1)*ncol + (jj-1)] = Ke_AB.yy(); Ke[(ii-1)*ncol + jj] = Ke_AB.yz();
          Ke[(ii-0)*ncol + (jj-2)] = Ke_AB.zx(); Ke[(ii-0)*ncol + (jj-1)] = Ke_AB.zy(); Ke[(ii-0)*ncol + jj] = Ke_AB.zz();
        }

        Vector_3 Ge_A = Vector_3( Ge[ii-2], Ge[ii-1], Ge[ii] );
        Vector_3 rot_Ge_A = QT.VecMult( Ge_A );  // rot_Ge_A = QT * Ge_A

        // Update Ge
        Ge[ii-2] = rot_Ge_A.x(); Ge[ii-1] = rot_Ge_A.y(); Ge[ii] = rot_Ge_A.z();
      }
    } // end loop over rows
  }
  else
    SYS_T::print_fatal("Error: this ringbc_type is not supported in PGAssem_Tet_CMM_GenAlpha.\n");
}


void PGAssem_Tet_Wall::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &sol_wall_disp,
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementw,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_wall_part,
    const IGenBC * const &gbc )
{
  // Residual & tangent contributions from the thin-walled linear membrane in CMM
  WallMembrane_KG( curr_time, dt, sol_a, sol_b, sol_wall_disp, lassem_ptr, elementw, quad_s, nbc_part, ringnbc_part, ebc_wall_part );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );
  
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_Wall::WallMembrane_KG(
    const double &curr_time,
    const double &dt, 
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const PDNSolution * const &sol_wall_disp,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_w,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_wall_part )
{
  const int dof_disp = 3; 

  double * array_a    = new double [nlgn * dof_mat ];
  double * array_b    = new double [nlgn * dof_mat ];
  double * array_c    = new double [nlgn * dof_disp];
  double * local_as   = new double [snLocBas * dof_mat ];
  double * local_bs   = new double [snLocBas * dof_mat ];
  double * local_cs   = new double [snLocBas * dof_disp];
  int    * LSIEN      = new    int [snLocBas];
  double * sctrl_x    = new double [snLocBas];
  double * sctrl_y    = new double [snLocBas];
  double * sctrl_z    = new double [snLocBas];
  double * sthickness = new double [snLocBas];
  double * syoungsmod = new double [snLocBas];
  double * sspringconst  = new double [snLocBas];
  double * sdampingconst = new double [snLocBas];
  double * quaprestress  = new double [ 6 * quad_s->get_num_quadPts() ];
  PetscInt * srow_index  = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  sol_wall_disp->GetLocalArray( array_c );

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_sele = ebc_wall_part -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_wall_part -> get_SIEN(ebc_id, ee, LSIEN);
    ebc_wall_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
    ebc_wall_part -> get_thickness(ee, sthickness);
    ebc_wall_part -> get_youngsmod(ee, syoungsmod);
    ebc_wall_part -> get_springconst( ee, sspringconst );
    ebc_wall_part -> get_dampingconst(ee, sdampingconst);
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_a, LSIEN, snLocBas, local_as);
    GetLocal(array_b, LSIEN, snLocBas, local_bs);

    GetLocal(array_c, LSIEN, snLocBas, dof_disp, local_cs);

    lassem_ptr->Assem_Tangent_Residual_EBC_Wall( curr_time, dt, local_as, local_bs, local_cs,
        element_w, sctrl_x, sctrl_y, sctrl_z, sthickness, syoungsmod,
        sspringconst, sdampingconst, quaprestress, quad_s);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
    }

    RingBC_KG( ringnbc_part, dof_mat, dof_mat * snLocBas, dof_mat * snLocBas,
        srow_index, srow_index, lassem_ptr->sur_Tangent, lassem_ptr->sur_Residual );

    MatSetValues(K, dof_mat*snLocBas, srow_index, dof_mat*snLocBas, srow_index,
        lassem_ptr->sur_Tangent, ADD_VALUES);

    VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
  }

  delete [] array_a;  array_a  = nullptr;
  delete [] array_b;  array_b  = nullptr;
  delete [] array_c;  array_c  = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] local_cs; local_cs = nullptr;
  delete [] LSIEN;    LSIEN    = nullptr;
  delete [] sctrl_x;  sctrl_x  = nullptr;
  delete [] sctrl_y;  sctrl_y  = nullptr;
  delete [] sctrl_z;  sctrl_z  = nullptr;
  delete [] sthickness;    sthickness    = nullptr;
  delete [] syoungsmod;    syoungsmod    = nullptr;
  delete [] sspringconst;  sspringconst  = nullptr;
  delete [] sdampingconst; sdampingconst = nullptr;
  delete [] quaprestress;  quaprestress  = nullptr;
  delete [] srow_index;    srow_index    = nullptr;
}


void PGAssem_Tet_Wall::Update_Wall_Prestress(
    const PDNSolution * const &sol_wall_disp,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_w,
    const IQuadPts * const &quad_s,
    ALocal_EBC * const &ebc_wall_part )
{
  const int dof_disp = 3;
  const int face_nqp = quad_s -> get_num_quadPts();

  double * array_b      = new double [nlgn * dof_disp];
  double * local_bs     = new double [snLocBas * dof_disp];
  int    * LSIEN        = new    int [snLocBas];
  double * sctrl_x      = new double [snLocBas];
  double * sctrl_y      = new double [snLocBas];
  double * sctrl_z      = new double [snLocBas];
  double * syoungsmod   = new double [snLocBas];
  double * quaprestress = new double [6 * face_nqp];

  std::vector<Matrix_3x3> sigma; sigma.resize( face_nqp );

  sol_wall_disp->GetLocalArray( array_b );

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_sele = ebc_wall_part -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_wall_part -> get_SIEN(ebc_id, ee, LSIEN);
    ebc_wall_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
    ebc_wall_part -> get_youngsmod(ee, syoungsmod  );
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_b, LSIEN, snLocBas, dof_disp, local_bs);

    element_w -> buildBasis( quad_s, sctrl_x, sctrl_y, sctrl_z );

    lassem_ptr->get_Wall_CauchyStress( local_bs, element_w, syoungsmod, sigma );

    // update prestress in Voigt notation (comps 11, 22, 33, 23, 13, 12)
    for(int qua=0; qua<face_nqp; ++qua)
    {
      quaprestress[6*qua]   += sigma[qua].xx();
      quaprestress[6*qua+1] += sigma[qua].yy();
      quaprestress[6*qua+2] += sigma[qua].zz();
      quaprestress[6*qua+3] += sigma[qua].yz();
      quaprestress[6*qua+4] += sigma[qua].xz();
      quaprestress[6*qua+5] += sigma[qua].xy();
    }

    ebc_wall_part -> set_prestress(ee, quaprestress);
  }

  delete [] array_b;      array_b      = nullptr;
  delete [] local_bs;     local_bs     = nullptr;
  delete [] LSIEN;        LSIEN        = nullptr;
  delete [] syoungsmod;   syoungsmod   = nullptr;
  delete [] quaprestress; quaprestress = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
}

// EOF
