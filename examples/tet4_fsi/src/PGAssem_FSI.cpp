#include "PGAssem_FSI.hpp"

PGAssem_FSI::PGAssem_FSI( 
    IPLocAssem_2x2Block * const &locassem_f_ptr,
    IPLocAssem_2x2Block * const &locassem_s_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    const IAGlobal_Mesh_Info * const &agmi_v,
    const IAGlobal_Mesh_Info * const &agmi_p,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_v,
    const ALocal_IEN * const &aien_p,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const ALocal_NodalBC * const &part_nbc_v,
    const ALocal_NodalBC * const &part_nbc_p,
    const ALocal_EBC * const &part_ebc,
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( agmi_v -> get_nLocBas() ),
  snLocBas(3),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn_v( pnode_v -> get_nlocghonode() ),
  nlgn_p( pnode_p -> get_nlocghonode() )
{
  // Make sure the data structure is compatible
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id)
  {
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_FSI, snLocBas has to be uniform.\n");
  }

  const int nloc_v = pnode_v -> get_nlocalnode();
  const int nloc_p = pnode_p -> get_nlocalnode();

  const int nlocrow_v = 3 * nloc_v;
  const int nlocrow_p = 1 * nloc_p;
  const int nlocrow = nlocrow_v + nlocrow_p;

  SYS_T::commPrint("     Empirical nonzero estimate: %d \n", in_nz_estimate);
  
  // Create matrix with routh preallocation
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 4*in_nz_estimate, PETSC_NULL, 4*in_nz_estimate, PETSC_NULL, &K);

  // Create vector
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  
  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();
  
  //Assem_nonzero_estimate( alelem_ptr, locassem_ptr, aien_ptr, pnode_ptr, part_nbc );

}


PGAssem_FSI::~PGAssem_FSI()
{
  VecDestroy(&G);
  MatDestroy(&K);
}


void PGAssem_FSI::EssBC_KG( const ALocal_NodalBC * const &nbc_v,
  const ALocal_NodalBC * const &nbc_p )
{
  // For three velocity fields
  for(int field=0; field<3; ++field)
  {
    const int local_dir = nbc_v -> get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_v -> get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
        MatSetValue(K, row, row, 1.0, ADD_VALUES);
      }
    }
    
    const int local_sla = nbc_v -> get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_v -> get_LPSN(field, ii);
        const int col = nbc_v -> get_LPMN(field, ii);
        MatSetValue(K, row, col, 1.0, ADD_VALUES);
        MatSetValue(K, row, row, -1.0, ADD_VALUES);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }
  
  // For pressure field
  const int local_dir = nbc_p -> get_Num_LD(0);
  if(local_dir > 0)
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      int row = nbc_p -> get_LDN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  } 
  
  const int local_sla = nbc_p -> get_Num_LPS(0);
  if(local_sla > 0)
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      int row = nbc_p -> get_LPSN(0, ii);
      int col = nbc_p -> get_LPMN(0, ii);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_FSI::EssBC_G( const ALocal_NodalBC * const &nbc_v, 
    const ALocal_NodalBC * const &nbc_p )
{
  // For three velocity fields
  for(int field=0; field<3; ++field)
  {
    const int local_dir = nbc_v -> get_Num_LD(field);
    if( local_dir > 0 )
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        const int row = nbc_v -> get_LDN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }

    const int local_sla = nbc_v -> get_Num_LPS(field);
    if( local_sla > 0 )
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        const int row = nbc_v -> get_LPSN(field, ii);
        VecSetValue(G, row, 0.0, INSERT_VALUES);
      }
    }
  }

  // For pressure field
  const int local_dir = nbc_p -> get_Num_LD(0);
  if(local_dir > 0)
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      int row = nbc_p -> get_LDN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_p -> get_Num_LPS(0);
  if(local_sla > 0)
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      int row = nbc_p -> get_LPSN(0, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

double PGAssem_FSI::Assem_surface_flowrate(
    const PDNSolution * const &disp,
    const PDNSolution * const &velo,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array_v = new double [nlgn_v * 3];
  double * array_d = new double [nlgn_v * 3];
  double * local_v = new double [snLocBas * 3];
  double * local_d = new double [snLocBas * 3];

  disp -> GetLocalArray( array_d );
  velo -> GetLocalArray( array_v ); 

  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN );

    ebc_part -> get_ctrlPts_xyz( ebc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    GetLocal( array_d, LSIEN, snLocBas, 3, local_d );
    GetLocal( array_v, LSIEN, snLocBas, 3, local_v );  

    esum += lassem_ptr -> get_flowrate( local_d, local_v, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array_v; delete [] array_d; delete [] local_v; delete [] local_d;
  array_v = nullptr; array_d = nullptr; local_v = nullptr; local_d = nullptr;

  delete [] LSIEN; delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  LSIEN = nullptr; sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum = 0.0;

  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}

double PGAssem_FSI::Assem_surface_ave_pressure(
    const PDNSolution * const &disp,
    const PDNSolution * const &pres,
    IPLocAssem_2x2Block * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array_d = new double [nlgn_v * 3];
  double * array_p = new double [nlgn_v];
  double * local_d = new double [snLocBas * 3];
  double * local_p = new double [snLocBas];
  
  disp -> GetLocalArray( array_d );
  pres -> GetLocalArray( array_p );
  
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;
  
  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN );

    ebc_part -> get_ctrlPts_xyz( ebc_id, ee, sctrl_x, sctrl_y, sctrl_z );

    GetLocal( array_d, LSIEN, snLocBas, 3, local_d );
    GetLocal( array_p, LSIEN, snLocBas, 1, local_p );

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local_d, local_p, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);
    
    val_pres += ele_pres;
    val_area += ele_area; 
  }

  delete [] array_p; delete [] array_d; delete [] local_p; delete [] local_d;
  array_p = nullptr; array_d = nullptr; local_p = nullptr; local_d = nullptr;

  delete [] LSIEN; delete [] sctrl_x; delete [] sctrl_y; delete [] sctrl_z;
  LSIEN = nullptr; sctrl_x = nullptr; sctrl_y = nullptr; sctrl_z = nullptr;

  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}



























// EOF
