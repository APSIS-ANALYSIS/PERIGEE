#include "PostVectSolution.hpp"

PostVectSolution::PostVectSolution( const std::string &solution_file_name,
    const std::string &analysis_node_mapping_file,
    const std::string &post_node_mapping_file,
    const APart_Node * const &aNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof )
: dof_per_node( input_dof ), 
  loc_sol_size( aNode_ptr->get_nlocghonode() * dof_per_node )
{
  loc_solution = new double [loc_sol_size];

  int nFunc = gInfo_ptr->get_nFunc();
  int whole_vec_size = nFunc * dof_per_node;
  double * vec_temp = new double [whole_vec_size];

  int * analysis_old2new = new int [nFunc];
  int * postproc_new2old = new int [nFunc];

  // Read PETSc solution vector into vec_temp
  ReadPETSc_vec(solution_file_name, whole_vec_size, vec_temp);

  // Read new2old and old2new mappings
  ReadNodeMapping(analysis_node_mapping_file, "old_2_new", nFunc, analysis_old2new );
  ReadNodeMapping(post_node_mapping_file, "new_2_old", nFunc, postproc_new2old );

  int index;
  for( int ii=0; ii<aNode_ptr->get_nlocghonode(); ++ii )
  {
    index = aNode_ptr->get_local_to_global(ii); // in postprocess partition's new index
    index = postproc_new2old[index]; // map back to natural global index
    index = analysis_old2new[index]; // map forward to analysis partitioned new index

    for(int jj=0; jj<dof_per_node; ++jj)
      loc_solution[ii*dof_per_node + jj] = vec_temp[index*dof_per_node + jj];
  }

  delete [] analysis_old2new;
  delete [] postproc_new2old;
  delete [] vec_temp;
}


PostVectSolution::PostVectSolution( const std::string &solution_file_name,
    const std::string &analysis_node_mapping_file,
    const std::string &post_node_mapping_file,
    const APart_Node * const &aNode_ptr,
    const int &in_nfunc, const int &input_dof )
: dof_per_node( input_dof ), 
  loc_sol_size( aNode_ptr->get_nlocghonode() * dof_per_node )
{
  loc_solution = new double [loc_sol_size];

  int nFunc = in_nfunc;
  int whole_vec_size = nFunc * dof_per_node;
  double * vec_temp = new double [whole_vec_size];

  int * analysis_old2new = new int [nFunc];
  int * postproc_new2old = new int [nFunc];

  // Read PETSc solution vector into vec_temp
  ReadPETSc_vec(solution_file_name, whole_vec_size, vec_temp);

  // Read new2old and old2new mappings
  ReadNodeMapping(analysis_node_mapping_file, "old_2_new", nFunc, analysis_old2new );
  ReadNodeMapping(post_node_mapping_file, "new_2_old", nFunc, postproc_new2old );

  int index;
  for( int ii=0; ii<aNode_ptr->get_nlocghonode(); ++ii )
  {
    index = aNode_ptr->get_local_to_global(ii); // in postprocess partition's new index
    index = postproc_new2old[index]; // map back to natural global index
    index = analysis_old2new[index]; // map forward to analysis partitioned new index

    for(int jj=0; jj<dof_per_node; ++jj)
      loc_solution[ii*dof_per_node + jj] = vec_temp[index*dof_per_node + jj];
  }

  delete [] analysis_old2new;
  delete [] postproc_new2old;
  delete [] vec_temp;
}



PostVectSolution::~PostVectSolution()
{
  delete [] loc_solution;
}


void PostVectSolution::print_info() const
{
  std::cout<<"loc_sol_size: "<<loc_sol_size<<std::endl;
  std::cout<<"dof_per_node: "<<dof_per_node<<std::endl;
  for(int ii=0; ii<loc_sol_size; ++ii)
    std::cout<<ii<<'\t'<<loc_solution[ii]<<'\n';
  std::cout<<std::endl;
}

void PostVectSolution::PlusAX(const PostVectSolution * const &input_sol, 
    const double &val)
{
  if(input_sol->get_dof() != dof_per_node)
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: dof of the two PostVectSolution does not match. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  if(input_sol->get_solsize() != loc_sol_size)
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: size of the two PostVectSolution does not match. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  for(int ii=0; ii<loc_sol_size; ++ii)
    loc_solution[ii] = loc_solution[ii] + val * input_sol->get_locsol(ii);
}


void PostVectSolution::get_esol(const int &field, const int &nLocBas,
    const int * const &eien, double * const &esol) const
{
  // check the input field index
  if(field >= dof_per_node)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Error: field is out of range. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  int index;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    index = eien[ii];
    esol[ii] = loc_solution[dof_per_node * index + field];
  }
}


void PostVectSolution::ReadPETSc_vec( const std::string &solution_file_name,
    const int &vec_size, double * const &veccopy )
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(),
      FILE_MODE_READ, &viewer);
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
  hid_t file_id = H5Fopen(node_mapping_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT); 

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0]; 

  if( int(dSize) != node_size )
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, nodemap );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);
}

// EOF
