#include "VisDataPrep_Block_HED_upv_3D.hpp"

VisDataPrep_Block_HED_upv_3D::VisDataPrep_Block_HED_upv_3D( const bool &is_ref )
{
  if(is_ref)
  {
    arrayCompSize = 3;

    arrayNames.push_back("Displacement");
    arraySizes.push_back(3);

    arrayNames.push_back("Pressure");
    arraySizes.push_back(1);

    arrayNames.push_back("Velocity");
    arraySizes.push_back(3);
  }
  else
  {
    arrayCompSize = 4;

    arrayNames.push_back("Displacement");
    arraySizes.push_back(3);
    
    arrayNames.push_back("detF");
    arraySizes.push_back(1);
    
    arrayNames.push_back("Pressure");
    arraySizes.push_back(1);

    arrayNames.push_back("Velocity");
    arraySizes.push_back(3);
  }

  pt_array_len.clear();
  pt_array_len.push_back(3);
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}



VisDataPrep_Block_HED_upv_3D::~VisDataPrep_Block_HED_upv_3D()
{}



void VisDataPrep_Block_HED_upv_3D::get_pointArray(
    const std::vector<std::string> solution_file_names,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const int &in_nfunc,
    double ** &pointArrays ) const
{
  int ntotal = nNode_ptr -> get_nlocghonode();

  PostVectSolution pvdisp( solution_file_names[0], 
      analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, in_nfunc, 3 );

  PostVectSolution pvpres( solution_file_names[1], 
      analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, in_nfunc, 1 );

  PostVectSolution pvvelo( solution_file_names[2], 
      analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, in_nfunc, 3 );

  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][3*ii]   = pvdisp.get_locsol(ii*3);
    pointArrays[0][3*ii+1] = pvdisp.get_locsol(ii*3+1);
    pointArrays[0][3*ii+2] = pvdisp.get_locsol(ii*3+2);
    pointArrays[1][ii]     = pvpres.get_locsol(ii);
    pointArrays[2][3*ii]   = pvvelo.get_locsol(ii*3);
    pointArrays[2][3*ii+1] = pvvelo.get_locsol(ii*3+1);
    pointArrays[2][3*ii+2] = pvvelo.get_locsol(ii*3+2);
  }

  if( get_ptarray_size() != 3) SYS_T::print_fatal("Error: get_ptarray_size != 3\n");
}


// EOF
