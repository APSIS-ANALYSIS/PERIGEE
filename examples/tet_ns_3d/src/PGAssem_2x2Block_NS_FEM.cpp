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

  // Temporarily ignore new entry allocation
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  // Allocate containers
  row_index_v = new PetscInt [nLocBas * dof_mat_v];
  row_index_p = new PetscInt [nLocBas * dof_mat_p];

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
    srow_index_v = new PetscInt [snLocBas * dof_mat_v];
    srow_index_p = new PetscInt [snLocBas * dof_mat_p];
  }

  // Initialize values for array
  for(int ii=0; ii<nlgn*dof_sol; ++ii)
  {
    array_a[ii] = 0.0;
    array_b[ii] = 0.0;
  }

  // Nonzero pattern assembly
  //Assem_nonzero_estimate( alelem_ptr, locassem_ptr,
  //    elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ebc, gbc );

}


PGAssem_2x2Block_NS_FEM::~PGAssem_2x2Block_NS_FEM()
{}


void PGAssem_2x2Block_NS_FEM::EssBC_KG( const ALocal_NodalBC * const &nbc_part )
{
  // pressure dof comes from field 0, to be inserted in subK[3] and subG[1]
  const int local_dir = nbc_part->get_Num_LD(0);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc_part->get_LDN(0, i) * dof_mat_p;
      MatSetValue(subK[0], row, row, 1.0, ADD_VALUES);
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(0);
  if(local_sla > 0)
  {
    for(int i=0; i<local_sla; ++i)
    {
      const int row = nbc_part->get_LPSN(0, i) * dof_mat_p;
      const int col = nbc_part->get_LPMN(0, i) * dof_mat_p;
      MatSetValue(subK[0], row, col, 1.0, ADD_VALUES);
      MatSetValue(subK[0], row, row, -1.0, ADD_VALUES);
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  // velocity dofs
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc_part->get_Num_LD(field);

    if(local_dir > 0)
    {
      for(int i=0; i<local_dir; ++i)
      {
        const int row = nbc_part->get_LDN(field, i) * dof_mat_v + field - 1;
        MatSetValue(subK[3], row, row, 1.0, ADD_VALUES);
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc_part->get_Num_LPS(field);
    if(local_sla > 0)
    {
      for(int i=0; i<local_sla; ++i)
      {
        const int row = nbc_part->get_LPSN(field, i) * dof_mat_v + field - 1;
        const int col = nbc_part->get_LPMN(field, i) * dof_mat_v + field - 1;
        MatSetValue(subK[3], row, col, 1.0, ADD_VALUES);
        MatSetValue(subK[3], row, row, -1.0, ADD_VALUES);
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
      }
    }
  }
}


void PGAssem_2x2Block_NS_FEM::EssBC_G( const ALocal_NodalBC * const &nbc_part )
{
  // pres field is 0, to be inserted to subG[1]
  const int local_dir = nbc_part->get_Num_LD(0);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc_part->get_LDN(0, ii) * dof_mat_p;
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(0);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc_part->get_LPSN(0, ii) * dof_mat_p;
      VecSetValue(subG[0], row, 0.0, INSERT_VALUES);
    }
  }

  // velo fields from 1 to 3
  for(int field=1; field<=3; ++field)
  {
    const int local_dir = nbc_part->get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_part->get_LDN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc_part->get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_part->get_LPSN(field, ii) * dof_mat_v + field - 1;
        VecSetValue(subG[1], row, 0.0, INSERT_VALUES);
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




void PGAssem_2x2Block_NS_FEM::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
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

      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat_v; ++mm)
          srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc_part -> get_LID(mm+1, LSIEN[ii]) + mm;
      }

      VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, lassem_ptr->Residual1, ADD_VALUES);
    }
  }
}


void PGAssem_2x2Block_NS_FEM::BackFlow_G( IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
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

      lassem_ptr->Assem_Residual_BackFlowStab( local_as, local_bs,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat_v; ++mm)
          srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc_part -> get_LID(mm+1, LSIEN[ii]) + mm;
      }

      VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, lassem_ptr->sur_Residual1, ADD_VALUES);
    }
  }
}


void PGAssem_2x2Block_NS_FEM::BackFlow_KG( const double &dt,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
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

      lassem_ptr->Assem_Tangent_Residual_BackFlowStab( dt, local_as, local_bs,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat_v; ++mm)
          srow_index_v[dof_mat_v * ii + mm] = dof_mat_v * nbc_part -> get_LID(mm+1, LSIEN[ii]) + mm;
      }

      MatSetValues(subK[3], dof_mat_v*snLocBas, srow_index_v, dof_mat_v*snLocBas, srow_index_v,
          lassem_ptr->sur_Tangent11, ADD_VALUES);

      VecSetValues(subG[1], dof_mat_v*snLocBas, srow_index_v, lassem_ptr->sur_Residual1, ADD_VALUES);
    }
  }
}


double PGAssem_2x2Block_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [vec -> get_nlgn()];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}


double PGAssem_2x2Block_NS_FEM::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_Inflow_NodalBC * const &infbc_part )
{
  double * array = new double [vec -> get_nlgn()];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell();

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( ee, LSIEN);

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}


double PGAssem_2x2Block_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [vec -> get_nlgn()];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}


double PGAssem_2x2Block_NS_FEM::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_Inflow_NodalBC * const &infbc_part )
{
  double * array = new double [ vec->get_nlgn() ];
  double * local = new double [snLocBas * dof_sol];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell();

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( ee, LSIEN);

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}


void PGAssem_2x2Block_NS_FEM::NatBC_Resis_G(
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id from solution vector dot_sol
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id from solution vector sol
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate );

    // P_n+alpha_f
    // lassem_ptr->get_model_para_1() gives alpha_f 
    const double val = P_n + lassem_ptr->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // Here, val is Pressure, and is used as the surface traction h = P I 
      // to calculate the boundary integral
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, val,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_index_v[3*ii+0] = dof_mat_v * nbc_part->get_LID(1, LSIEN[ii]);
        srow_index_v[3*ii+1] = dof_mat_v * nbc_part->get_LID(2, LSIEN[ii]) + 1;
        srow_index_v[3*ii+2] = dof_mat_v * nbc_part->get_LID(3, LSIEN[ii]) + 2;
      }

      VecSetValues(subG[1], snLocBas*3, srow_index_v, lassem_ptr->Residual1, ADD_VALUES);
    }
  }
}


void PGAssem_2x2Block_NS_FEM::NatBC_Resis_KG(
    const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const double a_f = lassem_ptr -> get_model_para_1();

  // dd_dv = dt x alpha_f x gamma
  const double dd_dv = dt * a_f * lassem_ptr->get_model_para_2();

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
    // Calculate dot flow rate for face with ebc_id and MPI_Allreduce them
    // Here, dot_sol is the solution at time step n+1 (not n+alpha_f!)
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id and MPI_Allreduce them
    // Here, sol is the solution at time step n+1 (not n+alpha_f!)
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate );

    // P_n+alpha_f 
    const double resis_val = P_n + a_f * (P_np1 - P_n);

    // Get the (potentially approximated) m := dP/dQ
    const double m_val = gbc -> get_m( ebc_id, dot_flrate, flrate );

    // Get the (potentially approximated) n := dP/d(dot_Q)
    const double n_val = gbc -> get_n( ebc_id, dot_flrate, flrate );
    
    // Define alpha_f * n + alpha_f * gamma * dt * m
    // coef a^t a enters as the consistent tangent for the resistance-type bc
    const double coef = a_f * n_val + dd_dv * m_val;

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
      Tan = nullptr;
      scol_idx = nullptr;
    }

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // For here, we scale the int_NA nx/y/z by factor 1
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, 1.0,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      // Residual vector is scaled by the resistance value
      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = resis_val * lassem_ptr->Residual1[3*ii];
        Res[3*ii+1] = resis_val * lassem_ptr->Residual1[3*ii+1];
        Res[3*ii+2] = resis_val * lassem_ptr->Residual1[3*ii+2];
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            // Residual1[3*A+ii] is intNB[A]*out_n[ii]
            Tan[temp_row + 3*B + 0] = coef * lassem_ptr->Residual1[3*A+ii] * intNB[B] * out_nx;
            Tan[temp_row + 3*B + 1] = coef * lassem_ptr->Residual1[3*A+ii] * intNB[B] * out_ny;
            Tan[temp_row + 3*B + 2] = coef * lassem_ptr->Residual1[3*A+ii] * intNB[B] * out_nz;
          }
        }
      }

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_idx[3*ii+0] = dof_mat_v * nbc_part->get_LID(1,LSIEN[ii]);
        srow_idx[3*ii+1] = dof_mat_v * nbc_part->get_LID(2,LSIEN[ii]) + 1;
        srow_idx[3*ii+2] = dof_mat_v * nbc_part->get_LID(3,LSIEN[ii]) + 2;
      }

      for(int ii=0; ii<num_face_nodes; ++ii)
      {
        scol_idx[ii*3+0] = dof_mat_v * map_Bj[ii*3+0];
        scol_idx[ii*3+1] = dof_mat_v * map_Bj[ii*3+1] + 1;
        scol_idx[ii*3+2] = dof_mat_v * map_Bj[ii*3+2] + 2;
      }

      MatSetValues(subK[3], snLocBas*3, srow_idx, num_face_nodes*3, scol_idx, Tan, ADD_VALUES);
      VecSetValues(subG[1], snLocBas*3, srow_idx, Res, ADD_VALUES);
    }

    if( num_face_nodes > 0 ) 
    {
      delete [] Tan; Tan = nullptr;
      delete [] scol_idx; scol_idx = nullptr;
    }
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
}


// EOF
