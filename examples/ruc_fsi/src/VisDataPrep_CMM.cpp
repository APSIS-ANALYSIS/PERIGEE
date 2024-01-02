#include "VisDataPrep_CMM.hpp"

VisDataPrep_CMM::VisDataPrep_CMM()
{
  // Data to be written
  arrayCompSize = 3;

  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
  arrayNames.push_back("WallDisp");
  arraySizes.push_back(3);

  // Data to be read
  pt_array_len.clear();
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
  pt_array_len.push_back(3);
}

void VisDataPrep_CMM::get_pointArray(
    const std::vector<std::string> solution_file_names,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const int &in_nfunc,
    double ** &solArrays ) const
{
  // Check that solution_file_names contains both pres/velo & disp inputs 
  if(solution_file_names.size() != 2)
    SYS_T::print_fatal("VisDataPrep_CMM Error: solution_file_names.size() != 2. \n");

  // pres, velo
  PostVectSolution pvsolu(solution_file_names[0], analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, in_nfunc, 4);

  // wall disp
  PostVectSolution pvdisp(solution_file_names[1], analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, in_nfunc, 3);

  // Total number of nodes to be read from the solution vector
  const int ntotal = nNode_ptr->get_nlocghonode();
 
  // Assign solution values to the corresponding physical field container
  for(int ii=0; ii<ntotal; ++ii)
  {
    solArrays[0][ii]     = pvsolu.get_locsol(ii*4+0);
    solArrays[1][3*ii]   = pvsolu.get_locsol(ii*4+1);
    solArrays[1][3*ii+1] = pvsolu.get_locsol(ii*4+2);
    solArrays[1][3*ii+2] = pvsolu.get_locsol(ii*4+3);
    solArrays[2][3*ii]   = pvdisp.get_locsol(ii*3+0);
    solArrays[2][3*ii+1] = pvdisp.get_locsol(ii*3+1);
    solArrays[2][3*ii+2] = pvdisp.get_locsol(ii*3+2);
  }

  // Check to make sure that ptarray_size gives correct output  
  if(get_ptarray_size() != 3)
    SYS_T::print_fatal("VisDataPrep_CMM Error: get_ptarray_size != 3. \n");
}

// EOF
