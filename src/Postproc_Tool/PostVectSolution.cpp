#include "PostVectSolution.hpp"

PostVectSolution::PostVectSolution( const std::string &solution_file_name,
    const std::vector<int> &analysis_node_mapping,
    const std::vector<int> &post_node_mapping,
    const APart_Node * const &aNode_ptr,
    int nFunc, int input_dof )
: dof_per_node( input_dof ), 
  loc_sol_size( aNode_ptr->get_nlocghonode() * dof_per_node ),
  loc_solution(loc_sol_size, 0.0)
{
  // Read the full PETSc solution vector into vec_temp
  const auto vec_temp = VIS_T::readPETSc_vec(solution_file_name);

  for( int ii=0; ii<aNode_ptr->get_nlocghonode(); ++ii )
  {
    int index = aNode_ptr->get_local_to_global(ii); // in postprocess partition's new index
    index = post_node_mapping[index];               // map back to natural global index
    index = analysis_node_mapping[index];           // map forward to analysis partitioned new index

    for(int jj=0; jj<dof_per_node; ++jj)
      loc_solution[ii*dof_per_node + jj] = vec_temp[index*dof_per_node + jj];
  }
}

void PostVectSolution::print_info() const
{
  std::cout<<"loc_sol_size: "<<loc_sol_size<<std::endl;
  std::cout<<"dof_per_node: "<<dof_per_node<<std::endl;
  for(int ii=0; ii<loc_sol_size; ++ii)
    std::cout<<ii<<'\t'<<loc_solution[ii]<<'\n';
  std::cout<<std::endl;
}

void PostVectSolution::get_esol( int field, int nLocBas,
    const int * const &eien, double * const &esol) const
{
  // check the input field index
  SYS_T::print_fatal_if( field >= dof_per_node, "Error: field is out of range. \n" );

  for(int ii=0; ii<nLocBas; ++ii)
    esol[ii] = loc_solution[dof_per_node * eien[ii] + field];
}

// EOF
