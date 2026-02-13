#include "VisDataPrep_NS.hpp"

VisDataPrep_NS::VisDataPrep_NS()
{
  // Data to be written
  arrayCompSize = 2;

  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);

  // Data to be read
  pt_array_len.clear();
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}

void VisDataPrep_NS::get_pointArray(
    const std::string solution_file_name,
    const std::vector<int> &analysis_node_mapping,
    const std::vector<int> &post_node_mapping,
    const APart_Node * const &nNode_ptr,
    const int &input_dof,
    double ** &solArrays ) const
{
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping,
      post_node_mapping, nNode_ptr, input_dof);

  // Total number of nodes to be read from the solution vector
  const int ntotal = nNode_ptr->get_nlocghonode();
 
  // Assign the solution values to the corresponding physical field
  // container 
  for(int ii=0; ii<ntotal; ++ii)
  {
    solArrays[0][ii]     = pvsolu.get_locsol(ii*input_dof+0);
    solArrays[1][3*ii]   = pvsolu.get_locsol(ii*input_dof+1);
    solArrays[1][3*ii+1] = pvsolu.get_locsol(ii*input_dof+2);
    solArrays[1][3*ii+2] = pvsolu.get_locsol(ii*input_dof+3);
  }

  // Check to make sure that ptarray_size gives correct output  
  if(get_ptarray_size() != 2) SYS_T::print_fatal("Error: get_ptarray_size != 2. \n");
}

// EOF
