#include "VisDataPrep_NS.hpp"
#include "PostVectSolution.hpp"
#include "APart_Node.hpp"

VisDataPrep_NS::VisDataPrep_NS(const bool &in_is_transport)
: is_transport(in_is_transport)
{
  // Data to be written
  arrayCompSize = is_transport ? 3 : 2;

  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
  if(is_transport)
  {
    arrayNames.push_back("HI");
    arraySizes.push_back(1);
  }

  // Data to be read
  pt_array_len.clear();
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
  if(is_transport) pt_array_len.push_back(1);
}

void VisDataPrep_NS::get_pointArray(
    const std::string solution_file_name,
    const std::vector<int> &analysis_node_mapping,
    const std::vector<int> &post_node_mapping,
    const APart_Node * const &nNode_ptr,
    double ** &solArrays ) const
{
  constexpr int dof = 4;

  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping,
      post_node_mapping, nNode_ptr, dof);

  // Total number of nodes to be read from the solution vector
  const int ntotal = nNode_ptr->get_nlocghonode();
 
  // Assign the solution values to the corresponding physical field
  // container 
  for(int ii=0; ii<ntotal; ++ii)
  {
    solArrays[0][ii]     = pvsolu.get_locsol(ii*dof+0);
    solArrays[1][3*ii]   = pvsolu.get_locsol(ii*dof+1);
    solArrays[1][3*ii+1] = pvsolu.get_locsol(ii*dof+2);
    solArrays[1][3*ii+2] = pvsolu.get_locsol(ii*dof+3);
  }

  // Check to make sure that ptarray_size gives correct output  
  if(get_ptarray_size() != 2) SYS_T::print_fatal("Error: get_ptarray_size != 2. \n");
}

void VisDataPrep_NS::get_pointArray(
    const std::vector<std::string> solution_file_names,
    const std::vector<int> &analysis_node_mapping,
    const std::vector<int> &post_node_mapping,
    const APart_Node * const &nNode_ptr,
    double ** &solArrays ) const
{
  SYS_T::print_fatal_if(solution_file_names.size() != 2,
      "Error: VisDataPrep_NS expects two solution files in hemolysis mode.\n");
  SYS_T::print_fatal_if(!is_transport,
      "Error: VisDataPrep_NS vector<string> reader requires hemolysis mode.\n");

  constexpr int flow_dof = 4;
  constexpr int scalar_dof = 1;

  PostVectSolution flow_sol(solution_file_names[0], analysis_node_mapping,
      post_node_mapping, nNode_ptr, flow_dof);
  PostVectSolution hi_sol(solution_file_names[1], analysis_node_mapping,
      post_node_mapping, nNode_ptr, scalar_dof);

  const int ntotal = nNode_ptr->get_nlocghonode();

  for(int ii=0; ii<ntotal; ++ii)
  {
    solArrays[0][ii] = flow_sol.get_locsol(ii*flow_dof+0);
    solArrays[1][3*ii] = flow_sol.get_locsol(ii*flow_dof+1);
    solArrays[1][3*ii+1] = flow_sol.get_locsol(ii*flow_dof+2);
    solArrays[1][3*ii+2] = flow_sol.get_locsol(ii*flow_dof+3);
    solArrays[2][ii] = hi_sol.get_locsol(ii*scalar_dof);
  }

  if(get_ptarray_size() != 3) SYS_T::print_fatal("Error: get_ptarray_size != 3. \n");
}

// EOF
