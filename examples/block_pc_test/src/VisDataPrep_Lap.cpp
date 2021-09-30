#include "VisDataPrep_Lap.hpp"

VisDataPrep_Lap::VisDataPrep_Lap()
{
  arrayCompSize = 1;

  arrayNames.push_back("u");
  arraySizes.push_back(1);

  pt_array_len.clear();
  pt_array_len.push_back(1);
}


VisDataPrep_Lap::~VisDataPrep_Lap()
{}


void VisDataPrep_Lap::get_pointArray(
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

  const int ntotal = nNode_ptr->get_nlocghonode();

  for(int ii=0; ii<ntotal; ++ii)
    solArrays[0][ii]     = pvsolu.get_locsol(ii*input_dof+0);

  if(get_ptarray_size() != 1) SYS_T::print_fatal("Error: get_ptarray_size != 1. \n");
}

// EOF
