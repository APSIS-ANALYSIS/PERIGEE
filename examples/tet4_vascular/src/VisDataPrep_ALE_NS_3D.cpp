#include "VisDataPrep_ALE_NS_3D.hpp"

VisDataPrep_ALE_NS_3D::VisDataPrep_ALE_NS_3D()
{
  arrayCompSize = 3;

  arrayNames.push_back("Displacement");
  arraySizes.push_back(3);
  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);

  pt_array_len.clear();
  pt_array_len.push_back(3);
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}


VisDataPrep_ALE_NS_3D::~VisDataPrep_ALE_NS_3D()
{}


void VisDataPrep_ALE_NS_3D::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  const int ntotal = nNode_ptr->get_nlocghonode();
  
  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][3*ii]   = pvsolu.get_locsol(ii*input_dof+0);
    pointArrays[0][3*ii+1] = pvsolu.get_locsol(ii*input_dof+1);
    pointArrays[0][3*ii+2] = pvsolu.get_locsol(ii*input_dof+2);
    pointArrays[1][ii]     = pvsolu.get_locsol(ii*input_dof+3);
    pointArrays[2][3*ii]   = pvsolu.get_locsol(ii*input_dof+4);
    pointArrays[2][3*ii+1] = pvsolu.get_locsol(ii*input_dof+5);
    pointArrays[2][3*ii+2] = pvsolu.get_locsol(ii*input_dof+6);
  }

  // Check to make sure that ptarray_size gives correct output  
  if(get_ptarray_size() != 3)
    SYS_T::print_fatal("Error: get_ptarray_size != 3. \n");
}

// EOF
