#include "Matrix_PETSc_CMM.hpp"

Matrix_PETSc_CMM::Matrix_PETSc_CMM( const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const ALocal_Ring_NodalBC * const &ring_bc_part )
: Matrix_PETSc( pnode_ptr, bc_part, 2, 0 )
{
  switch( ring_bc_part -> get_ringbc_type() )
  {
    case 0:
      gen_perm_bc( pnode_ptr, bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for clamped rings.\n");
      break;
    case 1:
      gen_ring_inplane_bc( pnode_ptr, bc_part, ring_bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for in-plane motion.\n");
      break;
    case 2:
      gen_ring_radial_motion_bc( pnode_ptr, bc_part, ring_bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for radial motion.\n");
      break;
    case 3:
      gen_outlet_ring_inplane_bc( pnode_ptr, bc_part, ring_bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for inlet clamping & outlet in-plane motion.\n");
      break;
    case 4:
      gen_outlet_ring_radial_motion_bc( pnode_ptr, bc_part, ring_bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for inlet clamping & outlet radial motion.\n");
      break;
    case 5:
      gen_ring_inplane_bc_partial_clamp( pnode_ptr, bc_part, ring_bc_part );
      SYS_T::commPrint("     Matrix_PETSc_CMM : matrix for in-plane motion & partial clamping.\n");
      break;
    default:
      SYS_T::print_fatal("Error: Matrix_PETSc_CMM has no such matrix type.\n");
      break;
  }
}


Matrix_PETSc_CMM::~Matrix_PETSc_CMM()
{}


void Matrix_PETSc_CMM::gen_ring_inplane_bc( const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_ring_inplane_bc assumes 4 dofs per node.\n");

  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      // For non-essential bc nodes, its diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();

  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
    const int dcomp = ring_bc_part -> get_dominant_n_comp( ii );
    const int row = dnode * dof + dcomp + 1;
    const double bot = -1.0 * ring_bc_part -> get_outvec(ii, dcomp);
    for(int jj=0; jj<3; ++jj)
    {
      if( jj != dcomp )
      {
        const int col = dnode * dof + jj + 1;
        const double top = ring_bc_part -> get_outvec(ii, jj);
        MatSetValue(K, row, col, top/bot, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc_CMM::gen_ring_radial_motion_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_ring_radial_motion_bc assumes 4 dofs per node.\n");

  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      // For non-essential bc nodes, its diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();

  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
    const int dncomp = ring_bc_part -> get_dominant_n_comp( ii );
    const int dtcomp = ring_bc_part -> get_dominant_t_comp( ii );
    const int nrow = dnode * dof + dncomp + 1;
    const int trow = dnode * dof + dtcomp + 1;

    // 3 - dncomp - dtcomp gives the remaining dof index. 
    const int rrow = dnode * dof + 4 - dncomp - dtcomp;

    const double na = ring_bc_part -> get_outvec(ii, dncomp);
    const double nb = ring_bc_part -> get_outvec(ii, dtcomp);
    const double nc = ring_bc_part -> get_outvec(ii, 3 - dncomp - dtcomp);
    
    const double ta = ring_bc_part -> get_tanvec(ii, dncomp);
    const double tb = ring_bc_part -> get_tanvec(ii, dtcomp);
    const double tc = ring_bc_part -> get_tanvec(ii, 3 - dncomp - dtcomp);
    
    const double val1 = (nb * tc - nc * tb) / (na * tb - nb * ta); 
    const double val2 = (ta * nc - na * tc) / (na * tb - nb * ta); 
        
    MatSetValue(K, nrow, rrow, val1, INSERT_VALUES);
    MatSetValue(K, trow, rrow, val2, INSERT_VALUES);
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc_CMM::gen_outlet_ring_inplane_bc( const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_outlet_ring_inplane_bc assumes 4 dofs per node.\n");

  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      // For non-essential bc nodes, its diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();

  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int cap_id = ring_bc_part -> get_cap_id( ii );

    if( cap_id != 0 )
    {
      const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
      const int dcomp = ring_bc_part -> get_dominant_n_comp( ii );
      const int row = dnode * dof + dcomp + 1;
      const double bot = -1.0 * ring_bc_part -> get_outvec(ii, dcomp);
      for(int jj=0; jj<3; ++jj)
      {
        if( jj != dcomp )
        {
          const int col = dnode * dof + jj + 1;
          const double top = ring_bc_part -> get_outvec(ii, jj);
          MatSetValue(K, row, col, top/bot, INSERT_VALUES);
        }
      }
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc_CMM::gen_outlet_ring_radial_motion_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_outlet_ring_radial_motion_bc assumes 4 dofs per node.\n");

  for(int ii=0; ii<nnode; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;
      // For non-essential bc nodes, its diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();

  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int cap_id = ring_bc_part -> get_cap_id( ii );

    if( cap_id != 0 )
    {
      const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
      const int dncomp = ring_bc_part -> get_dominant_n_comp( ii );
      const int dtcomp = ring_bc_part -> get_dominant_t_comp( ii );
      const int nrow = dnode * dof + dncomp + 1;
      const int trow = dnode * dof + dtcomp + 1;

      // 3 - dncomp - dtcomp gives the remaining dof index. 
      const int rrow = dnode * dof + 4 - dncomp - dtcomp;

      const double na = ring_bc_part -> get_outvec(ii, dncomp);
      const double nb = ring_bc_part -> get_outvec(ii, dtcomp);
      const double nc = ring_bc_part -> get_outvec(ii, 3 - dncomp - dtcomp);
      
      const double ta = ring_bc_part -> get_tanvec(ii, dncomp);
      const double tb = ring_bc_part -> get_tanvec(ii, dtcomp);
      const double tc = ring_bc_part -> get_tanvec(ii, 3 - dncomp - dtcomp);
      
      const double val1 = (nb * tc - nc * tb) / (na * tb - nb * ta); 
      const double val2 = (ta * nc - na * tc) / (na * tb - nb * ta); 
          
      MatSetValue(K, nrow, rrow, val1, INSERT_VALUES);
      MatSetValue(K, trow, rrow, val2, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}

void Matrix_PETSc_CMM::gen_ring_inplane_bc_partial_clamp( const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const ALocal_Ring_NodalBC * const &ring_bc_part )
{
  if(is_set) Clear();
  SYS_T::print_fatal_if(m != n, "Error: This is not a square matrix. \n");

  const int nnode = pnode_ptr->get_nlocalnode();
  const int dof   = bc_part->get_dof_LID();

  // We enforce 4 dofs per node
  SYS_T::print_fatal_if(dof != 4, "Error: Matrix_PETSc::gen_ring_inplane_bc_partial_clamp assumes 4 dofs per node.\n");

  std::vector<int> clamped_nodes;
  for(int ii=0; ii<nnode; ++ii)
  {
    int num_ess_velo_dof = 0;
    for(int jj=0; jj<dof; ++jj)
    {
      const int row = pnode_ptr->get_node_loc(ii) * dof + jj;
      const int col = bc_part->get_LID(jj, ii) * dof + jj;

      if( jj >= 1 && col < 0 ) num_ess_velo_dof += 1;

      // For non-essential bc nodes, its diagonal entry will be assigned 1.0
      MatSetValue(K, row, col, 1.0, INSERT_VALUES);
    }

    if( num_ess_velo_dof == 3 ) clamped_nodes.push_back( pnode_ptr->get_node_loc(ii) );
  }

  const int num_ring_node = ring_bc_part -> get_Num_LD();

  for(int ii=0; ii<num_ring_node; ++ii)
  {
    const int dnode = ring_bc_part -> get_LDN( ii ); // The ring node's nodal index
    const int dcomp = ring_bc_part -> get_dominant_n_comp( ii );
    const int row = dnode * dof + dcomp + 1;
    const double bot = -1.0 * ring_bc_part -> get_outvec(ii, dcomp);

    if( !VEC_T::is_invec( clamped_nodes, dnode ) )
    {
      for(int jj=0; jj<3; ++jj)
      {
        if( jj != dcomp )
        {
          const int col = dnode * dof + jj + 1;
          const double top = ring_bc_part -> get_outvec(ii, jj);
          MatSetValue(K, row, col, top/bot, INSERT_VALUES);
        }
      }
    }
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  is_set = true;
}


// EOF
