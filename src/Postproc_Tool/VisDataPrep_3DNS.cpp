#include "VisDataPrep_3DNS.hpp"

VisDataPrep_3DNS::VisDataPrep_3DNS()
{
  arrayCompSize = 2;
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
}

VisDataPrep_3DNS::~VisDataPrep_3DNS()
{}

void VisDataPrep_3DNS::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  int ntotal = nNode_ptr->get_nlocghonode();

  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][3*ii]   = pvsolu.get_locsol(ii*input_dof);
    pointArrays[0][3*ii+1] = pvsolu.get_locsol(ii*input_dof+1);
    pointArrays[0][3*ii+2] = pvsolu.get_locsol(ii*input_dof+2);
    pointArrays[1][ii]     = pvsolu.get_locsol(ii*input_dof+3);
  }
}

// EOF
