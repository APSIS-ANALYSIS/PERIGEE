#include "PostVectSolution.hpp"

PostVectSolution::PostVectSolution( const std::string &solution_file_name,
    const std::string &analysis_node_mapping_file,
    const std::string &post_node_mapping_file,
    const APart_Node * const &aNode_ptr,
    const int &nFunc, const int &input_dof )
: dof_per_node( input_dof ), 
  loc_sol_size( aNode_ptr->get_nlocghonode() * dof_per_node ),
  loc_solution(loc_sol_size, 0.0)
{
  double * vec_temp = new double [ nFunc * dof_per_node ];

  int * analysis_old2new = new int [nFunc];
  int * postproc_new2old = new int [nFunc];

  // Read the full PETSc solution vector into vec_temp
  ReadPETSc_vec(solution_file_name, nFunc * dof_per_node, vec_temp);

  // Read new2old and old2new mappings from HDF5 files
  ReadNodeMapping(analysis_node_mapping_file, "old_2_new", nFunc, analysis_old2new );
  ReadNodeMapping(post_node_mapping_file, "new_2_old", nFunc, postproc_new2old );

  for( int ii=0; ii<aNode_ptr->get_nlocghonode(); ++ii )
  {
    int index = aNode_ptr->get_local_to_global(ii); // in postprocess partition's new index
    index = postproc_new2old[index];                // map back to natural global index
    index = analysis_old2new[index];                // map forward to analysis partitioned new index

    for(int jj=0; jj<dof_per_node; ++jj)
      loc_solution[ii*dof_per_node + jj] = vec_temp[index*dof_per_node + jj];
  }

  delete [] analysis_old2new; analysis_old2new = nullptr;
  delete [] postproc_new2old; postproc_new2old = nullptr;
  delete [] vec_temp;         vec_temp         = nullptr;
}

void PostVectSolution::print_info() const
{
  std::cout<<"loc_sol_size: "<<loc_sol_size<<std::endl;
  std::cout<<"dof_per_node: "<<dof_per_node<<std::endl;
  for(int ii=0; ii<loc_sol_size; ++ii)
    std::cout<<ii<<'\t'<<loc_solution[ii]<<'\n';
  std::cout<<std::endl;
}

void PostVectSolution::get_esol(const int &field, const int &nLocBas,
    const int * const &eien, double * const &esol) const
{
  // check the input field index
  SYS_T::print_fatal_if( field >= dof_per_node, "Error: field is out of range. \n" );

  for(int ii=0; ii<nLocBas; ++ii)
    esol[ii] = loc_solution[dof_per_node * eien[ii] + field];
}

void PostVectSolution::ReadPETSc_vec( const std::string &solution_file_name,
    const int &vec_size, double * const &veccopy )
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(), FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the sol_temp has correct size
  PetscInt get_sol_temp_size;
  VecGetSize(sol_temp, &get_sol_temp_size);
  if( get_sol_temp_size != vec_size )
  {
    PetscPrintf(PETSC_COMM_SELF,
        "The solution size %d is not compatible with the size %d given by partition file! \n",
        get_sol_temp_size, vec_size);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  // read in array
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii)
    veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);
}

void PostVectSolution::ReadNodeMapping( const std::string &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap ) const
{
  const std::vector<int> temp_nodemap = HDF5_T::read_intVector( node_mapping_file.c_str(), "/", mapping_type );

  SYS_T::print_fatal_if(int( temp_nodemap.size() ) != node_size, "Error: PostVectSolution the allocated array has wrong size! \n");

  for(int ii=0; ii<node_size; ++ii) nodemap[ii] = temp_nodemap[ii];
}

// EOF
