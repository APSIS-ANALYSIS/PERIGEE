#include "VisDataPrep_Hyperelastic.hpp"

VisDataPrep_Hyperelastic::VisDataPrep_Hyperelastic( const bool &is_ref )
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
  arrayCompSize = 6;

  arrayNames.push_back("Displacement");
  arraySizes.push_back(3);
  arrayNames.push_back("detF");
  arraySizes.push_back(1);
  arrayNames.push_back("Pressure");
  arraySizes.push_back(1);
  arrayNames.push_back("Velocity");
  arraySizes.push_back(3);
  arrayNames.push_back("VonMisesStress");
  arraySizes.push_back(1);
  arrayNames.push_back("CauchyStress");
  arraySizes.push_back(6);
  }

  pt_array_len.clear();
  pt_array_len.push_back(3);
  pt_array_len.push_back(1);
  pt_array_len.push_back(3);
}

VisDataPrep_Hyperelastic::~VisDataPrep_Hyperelastic()
{}

void VisDataPrep_Hyperelastic::get_pointArray(
    const std::string &disp_sol_file_name,
    const std::string &pres_sol_file_name,
    const std::string &velo_sol_file_name,
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
  // Read local disp vector
  PostVectSolution pvsolu_disp(disp_sol_file_name, an_v_mapping_file,
      pn_v_mapping_file, pNode_v, input_nfunc_v, 3);

  // Read local velo vector
  PostVectSolution pvsolu_velo(velo_sol_file_name, an_v_mapping_file,
      pn_v_mapping_file, pNode_v, input_nfunc_v, 3);

  // Read local pres vector
  PostVectSolution pvsolu_pres(pres_sol_file_name, an_p_mapping_file,
      pn_p_mapping_file, pNode_p, input_nfunc_p, 1);

  const int ntotal_v = pNode_v -> get_nlocghonode();
  const int ntotal_p = pNode_p -> get_nlocghonode();

  for(int ii=0; ii<ntotal_v; ++ii)
  {
    pointArrays[0][3*ii  ] = pvsolu_disp.get_locsol(ii*3  );
    pointArrays[0][3*ii+1] = pvsolu_disp.get_locsol(ii*3+1);
    pointArrays[0][3*ii+2] = pvsolu_disp.get_locsol(ii*3+2);

    pointArrays[2][3*ii  ] = pvsolu_velo.get_locsol(ii*3  );
    pointArrays[2][3*ii+1] = pvsolu_velo.get_locsol(ii*3+1);
    pointArrays[2][3*ii+2] = pvsolu_velo.get_locsol(ii*3+2);
  }

  for(int ii=0; ii<ntotal_p; ++ii)
    pointArrays[1][ii]     = pvsolu_pres.get_locsol(ii);
}

// EOF
