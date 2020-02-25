#include "VisDataPrep_NS_3D.hpp"

VisDataPrep_NS_3D::VisDataPrep_NS_3D()
{
  arrayCompSize = 2;

  arrayNames.push_back("Pressure / rho_0");
  arraySizes.push_back(1);

  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);

  pt_array_len.clear();
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}

VisDataPrep_NS_3D::~VisDataPrep_NS_3D()
{}

void VisDataPrep_NS_3D::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  // Make sure the input sol has two types data to be read
  if(get_ptarray_size() != 2) SYS_T::print_fatal("Error: get_ptarray_size != 2. \n");
  
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  const int ntotal = nNode_ptr->get_nlocghonode();
  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][ii]     = pvsolu.get_locsol(ii*input_dof+0);
    pointArrays[1][3*ii]   = pvsolu.get_locsol(ii*input_dof+1);
    pointArrays[1][3*ii+1] = pvsolu.get_locsol(ii*input_dof+2);
    pointArrays[1][3*ii+2] = pvsolu.get_locsol(ii*input_dof+3);
  }
}

// EOF
