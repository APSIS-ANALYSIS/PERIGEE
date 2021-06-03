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
  const int local_dir = nbc_part->get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc_part->get_LDN(field, i) * dof_mat + field;

      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }
}


void PGAssem_Tet_Wall::EssBC_G( const ALocal_NodalBC * const &nbc_part, 
    const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc_part->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
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
  WallMembrane_KG( curr_time, dt, sol_a, sol_b, sol_wall_disp, lassem_ptr, elementw, quad_s, nbc_part, ebc_wall_part );

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
  double * quaprestress = new double [ 6 * quad_s->get_num_quadPts() ];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

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
    ebc_wall_part -> get_thickness(ee, sthickness  );
    ebc_wall_part -> get_youngsmod(ee, syoungsmod  );
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_a, LSIEN, snLocBas, local_as);
    GetLocal(array_b, LSIEN, snLocBas, local_bs);

    GetLocal(array_c, LSIEN, snLocBas, dof_disp, local_cs);

    lassem_ptr->Assem_Tangent_Residual_EBC_Wall( curr_time, dt, local_as, local_bs, local_cs,
        element_w, sctrl_x, sctrl_y, sctrl_z, sthickness, syoungsmod, quaprestress, quad_s);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
    }

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
  delete [] sthickness; sthickness = nullptr;
  delete [] syoungsmod; syoungsmod = nullptr;
  delete [] quaprestress; quaprestress = nullptr;
  delete [] srow_index; srow_index = nullptr;
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

    lassem_ptr->get_Wall_CauchyStress( local_bs, element_w, syoungsmod, quad_s, sigma );

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
