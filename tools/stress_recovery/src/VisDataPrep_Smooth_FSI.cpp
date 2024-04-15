#include "VisDataPrep_Smooth_FSI.hpp"

VisDataPrep_Smooth_FSI::VisDataPrep_Smooth_FSI()
{
  // Data to be written
  arrayCompSize = 8;

  arrayNames.push_back("CauchyStress");
  arraySizes.push_back(9);
  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Von Mises");
  arraySizes.push_back(1);
  arrayNames.push_back("Grad Disp");
  arraySizes.push_back(9);
  arrayNames.push_back("F");
  arraySizes.push_back(9);
  arrayNames.push_back("Von Mises_nop");
  arraySizes.push_back(1);
  arrayNames.push_back("Strain");
  arraySizes.push_back(9);
  arrayNames.push_back("disp");
  arraySizes.push_back(3);

  // Data to be read
  pt_array_len.clear();
  pt_array_len.push_back(9);
  pt_array_len.push_back(1);
  pt_array_len.push_back(1);
  pt_array_len.push_back(9);
  pt_array_len.push_back(9);
  pt_array_len.push_back(9);
  pt_array_len.push_back(9);
  pt_array_len.push_back(3);
}

void VisDataPrep_Smooth_FSI::get_pointArray(
    const std::string &grad_sol_file_name,
    const std::string &disp_sol_file_name,
    const std::string &pres_sol_file_name,
    const std::string &an_v_mapping_file,
    const std::string &an_p_mapping_file,
    const std::string &pn_v_mapping_file,
    const std::string &pn_p_mapping_file,
    const APart_Node * const &pNode_v,
    const APart_Node * const &pNode_p,
    const int &input_nfunc_v,
    const int &input_nfunc_p,
    double ** &pointArrays ) const
{
  PostVectSolution pvsolu_grad(grad_sol_file_name, an_v_mapping_file,
      pn_v_mapping_file, pNode_v, input_nfunc_v, 9);

  PostVectSolution pvsolu_pres(pres_sol_file_name, an_p_mapping_file,
      pn_p_mapping_file, pNode_p, input_nfunc_p, 1);

  PostVectSolution pvsolu_disp(disp_sol_file_name, an_v_mapping_file,
      pn_v_mapping_file, pNode_v, input_nfunc_v, 3);

  // Total number of nodes to be read from the solution vector
  const int ntotal_v = pNode_v -> get_nlocghonode();
  const int ntotal_p = pNode_p -> get_nlocghonode();
 
  // Assign the solution values to the corresponding physical field
  // container 
  for(int ii=0; ii<ntotal_v; ++ii)
  {
    pointArrays[0][9*ii]   = pvsolu_grad.get_locsol(9*ii+0);
    pointArrays[0][9*ii+1] = pvsolu_grad.get_locsol(9*ii+1);
    pointArrays[0][9*ii+2] = pvsolu_grad.get_locsol(9*ii+2);
    pointArrays[0][9*ii+3] = pvsolu_grad.get_locsol(9*ii+3);
    pointArrays[0][9*ii+4] = pvsolu_grad.get_locsol(9*ii+4);
    pointArrays[0][9*ii+5] = pvsolu_grad.get_locsol(9*ii+5);
    pointArrays[0][9*ii+6] = pvsolu_grad.get_locsol(9*ii+6);
    pointArrays[0][9*ii+7] = pvsolu_grad.get_locsol(9*ii+7);
    pointArrays[0][9*ii+8] = pvsolu_grad.get_locsol(9*ii+8);   
  }

  for(int ii=0; ii<ntotal_p; ++ii)
    pointArrays[1][ii] = pvsolu_pres.get_locsol(ii);

  for(int ii=0; ii<ntotal_v; ++ii)
  {
    pointArrays[2][3*ii]   = pvsolu_disp.get_locsol(3*ii+0);
    pointArrays[2][3*ii+1] = pvsolu_disp.get_locsol(3*ii+1);
    pointArrays[2][3*ii+2] = pvsolu_disp.get_locsol(3*ii+2);    
  }
}

// EOF
