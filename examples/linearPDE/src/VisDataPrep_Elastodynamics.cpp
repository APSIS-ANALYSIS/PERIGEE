#include "VisDataPrep_Elastodynamics.hpp"

VisDataPrep_Elastodynamics::VisDataPrep_Elastodynamics()
{
  // Data to be written
  arrayCompSize = 1;

  arrayNames.push_back("Displacement");
  arraySizes.push_back(3);

  // Data to be read
  pt_array_len.clear();
  pt_array_len.push_back(3);
}

void VisDataPrep_Elastodynamics::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const int &input_nfunc,
    const int &input_dof,
    double ** &solArrays ) const
{
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, input_nfunc, input_dof);

  // Total number of nodes to be read from the solution vector
  const int ntotal = nNode_ptr->get_nlocghonode();
 
  // Assign the solution values to the corresponding physical field
  // container 
  for(int ii=0; ii<ntotal; ++ii)
  {
    solArrays[0][3*ii]   = pvsolu.get_locsol(ii*input_dof+0);
    solArrays[0][3*ii+1] = pvsolu.get_locsol(ii*input_dof+1);
    solArrays[0][3*ii+2] = pvsolu.get_locsol(ii*input_dof+2);
  }

  // Check to make sure that ptarray_size gives correct output  
  if(get_ptarray_size() != 1) SYS_T::print_fatal("Error: get_ptarray_size != 1. \n");
}

// EOF
