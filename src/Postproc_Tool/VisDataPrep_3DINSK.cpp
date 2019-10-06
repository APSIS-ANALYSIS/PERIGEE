#include "VisDataPrep_3DINSK.hpp"

VisDataPrep_3DINSK::VisDataPrep_3DINSK()
{
  arrayCompSize = 2;
  arrayNames.push_back("Density");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
}


VisDataPrep_3DINSK::~VisDataPrep_3DINSK()
{}


void VisDataPrep_3DINSK::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  int ntotal = nNode_ptr->get_nlocghonode();
  
  // Read in the postprocessor's local solution from petsc binary
  // this constructor automatically handle the different index method
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][ii]     = pvsolu.get_locsol(ii*input_dof);
    pointArrays[1][3*ii]   = pvsolu.get_locsol(ii*input_dof+1);
    pointArrays[1][3*ii+1] = pvsolu.get_locsol(ii*input_dof+2);
    pointArrays[1][3*ii+2] = pvsolu.get_locsol(ii*input_dof+3);
  }
}

//EOF
