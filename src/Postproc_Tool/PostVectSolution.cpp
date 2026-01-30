#include "PostVectSolution.hpp"

PostVectSolution::PostVectSolution( const std::string &solution_file_name,
    const std::vector<int> &analysis_node_mapping,
    const std::vector<int> &post_node_mapping,
    const APart_Node * const &aNode_ptr,
    const int &nFunc, const int &input_dof )
: dof_per_node( input_dof ), 
  loc_sol_size( aNode_ptr->get_nlocghonode() * dof_per_node )
{
  // Allocate the space for the data loc_solution
  loc_solution = new double [loc_sol_size];

  // Read the full PETSc solution vector into vec_temp
  std::vector<double> vec_temp = VIS_T::readPETSc_vec(solution_file_name, nFunc * dof_per_node);

  for( int ii=0; ii<aNode_ptr->get_nlocghonode(); ++ii )
  {
    int index = aNode_ptr->get_local_to_global(ii); // in postprocess partition's new index
    index = post_node_mapping[index];                // map back to natural global index
    index = analysis_node_mapping[index];                // map forward to analysis partitioned new index

    for(int jj=0; jj<dof_per_node; ++jj)
      loc_solution[ii*dof_per_node + jj] = vec_temp[index*dof_per_node + jj];
  }
}

PostVectSolution::~PostVectSolution()
{
  delete [] loc_solution; loc_solution = nullptr;
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
