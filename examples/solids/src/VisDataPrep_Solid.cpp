#include "VisDataPrep_Solid.hpp"
#include "Sys_Tools.hpp"

VisDataPrep_Solid::VisDataPrep_Solid()
{
  arrayCompSize = 5;

  arrayNames.push_back("Displacement");
  arraySizes.push_back(3);
  arrayNames.push_back("detF");
  arraySizes.push_back(1);
  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
  arrayNames.push_back("CauchyStress");
  arraySizes.push_back(6);

  pt_array_len.clear();
  pt_array_len.push_back(3);
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}

void VisDataPrep_Solid::get_pointArray(
    const std::vector<std::string> solution_file_names,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const int &input_nfunc,
    double ** &pointArrays ) const
{
  SYS_T::print_fatal_if(solution_file_names.size() != 3,
      "Error: VisDataPrep_Solid expects three solution files (disp, pres, velo).\n");

  PostVectSolution pvsolu_disp(solution_file_names[0],
      analysis_node_mapping_file, post_node_mapping_file,
      nNode_ptr, input_nfunc, 3);

  PostVectSolution pvsolu_pres(solution_file_names[1],
      analysis_node_mapping_file, post_node_mapping_file,
      nNode_ptr, input_nfunc, 1);

  PostVectSolution pvsolu_velo(solution_file_names[2],
      analysis_node_mapping_file, post_node_mapping_file,
      nNode_ptr, input_nfunc, 3);

  const int ntotal = nNode_ptr->get_nlocghonode();

  for(int ii=0; ii<ntotal; ++ii)
  {
    pointArrays[0][3*ii  ] = pvsolu_disp.get_locsol(ii*3  );
    pointArrays[0][3*ii+1] = pvsolu_disp.get_locsol(ii*3+1);
    pointArrays[0][3*ii+2] = pvsolu_disp.get_locsol(ii*3+2);

    pointArrays[1][ii] = pvsolu_pres.get_locsol(ii);

    pointArrays[2][3*ii  ] = pvsolu_velo.get_locsol(ii*3  );
    pointArrays[2][3*ii+1] = pvsolu_velo.get_locsol(ii*3+1);
    pointArrays[2][3*ii+2] = pvsolu_velo.get_locsol(ii*3+2);
  }
}

// EOF
