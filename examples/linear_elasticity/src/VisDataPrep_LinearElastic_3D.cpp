#include "VisDataPrep_LinearElastic_3D.hpp"

VisDataPrep_LinearElastic_3D::VisDataPrep_LinearElastic_3D()
{
  arrayCompSize = 1;
  arrayNames.push_back("Displacement");
  arraySizes.push_back(3);

  pt_array_len.clear();
  pt_array_len.push_back(3);
}


VisDataPrep_LinearElastic_3D::~VisDataPrep_LinearElastic_3D()
{}


void VisDataPrep_LinearElastic_3D::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  int ntotal = nNode_ptr -> get_nlocghonode();
  
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
            post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][3*ii]   = pvsolu.get_locsol(ii*input_dof);
    pointArrays[0][3*ii+1] = pvsolu.get_locsol(ii*input_dof+1);
    pointArrays[0][3*ii+2] = pvsolu.get_locsol(ii*input_dof+2);
  }
}


// EOF
