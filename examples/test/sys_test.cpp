#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "ElemBC_3D_turbulence_wall_model.hpp"
#include "EBC_Partition_outflow.hpp"
#include "EBC_Partition_turbulence_wall_model.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include <iomanip>
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_vol_ST91.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "MaterialModel_NeoHookean_ST91_Mixed.hpp"
#include <memory>

int main(int argc, char *argv[])
{
  const double rho0 = 1.23;
  const double elastic_E = 1.2555; 
  const double elastic_nu = 0.47;
  IMaterialModel * oldmodel = new MaterialModel_NeoHookean_ST91_Mixed(rho0, elastic_E, elastic_nu);
  
  std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_ST91>(rho0, oldmodel->get_elastic_kappa());
  std::unique_ptr<IMaterialModel_ich> imodel = SYS_T::make_unique<MaterialModel_ich_NeoHookean>(oldmodel->get_elastic_mu());

  MaterialModel_Mixed_Elasticity * matmodel = new MaterialModel_Mixed_Elasticity(std::move(vmodel), std::move(imodel));

  matmodel->print_info();

  std::cout<<matmodel->get_model_name()<<'\n';

  //oldmodel->print_info();

  Tensor2_3D F; F.gen_rand();

  if(F.det() < 0.0) F*=-1.0;

  std::cout<<"F:"<<std::endl; 
  F.print();
  std::cout<<std::endl; 

  const double val = std::pow(F.det(), -1.0/3.0);

  //std::cout<<val<<'\n';

  F *= val; 
  
  Tensor2_3D P_old, S_old;
  oldmodel->get_PK(F, P_old, S_old);
  // S
  //S_old -= S_new;
  auto S_new = matmodel->get_PK_2nd(F);
  S_old -= S_new.full();
  S_old.print();
  std::cout<<std::endl;
  S_new.print();
  std::cout<<std::endl; 


  //P_old -= P_new
  auto P_new = matmodel->get_PK_1st(F);
  P_old -= P_new;
  P_old.print();
  std::cout<<std::endl;
  P_new.print();
  std::cout<<std::endl; 

  //CC_old -= CC_new
  Tensor4_3D CC_old;
  oldmodel->get_PK_Stiffness(F, P_old, S_old, CC_old);
  auto CC_new = matmodel->get_PK_Stiffness(F, P_new);
  CC_old -= CC_new.full();
  CC_old.print();
  std::cout<<std::endl;
  CC_new.print();
  std::cout<<std::endl; 

  double p = MATH_T::gen_double_rand(-1.0,1.0);
  // rho
  auto rho_old = oldmodel->get_rho(p);
  auto rho_new = matmodel->get_rho(p);
  rho_old -= rho_new;
  std::cout<<rho_old<<std::endl;
  std::cout<<rho_new<<std::endl;

  // beta
  auto beta_old = oldmodel->get_beta(p);
  auto beta_new = matmodel->get_beta(p);
  beta_old -= beta_new;
  std::cout<<beta_old<<std::endl;
  std::cout<<beta_new<<std::endl;

  // drho_dp
  auto drho_dp_old = oldmodel->get_drho_dp(p);
  auto drho_dp_new = matmodel->get_drho_dp(p);
  drho_dp_old -= drho_dp_new;
  std::cout<<drho_dp_old<<std::endl;
  std::cout<<drho_dp_new<<std::endl;

  // dbeta_dp
  auto dbeta_dp_old = oldmodel->get_dbeta_dp(p);
  auto dbeta_dp_new = matmodel->get_dbeta_dp(p);
  dbeta_dp_old -= dbeta_dp_new;
  std::cout<<dbeta_dp_old<<std::endl;
  std::cout<<dbeta_dp_new<<std::endl;

  // strain_energy
  auto strain_energy_old = oldmodel->get_strain_energy(F);
  auto strain_energy_Gibbs_new = matmodel->get_ich_energy(F) + matmodel->get_vol_Gibbs_energy(p);
  auto strain_energy_Helm_new = matmodel->get_ich_energy(F) + matmodel->get_vol_Helmholtz_energy(F.det());
  auto strain_energy_diff_Gibbs = strain_energy_old - strain_energy_Gibbs_new;
  auto strain_energy_diff_Helm = strain_energy_old - strain_energy_Helm_new;
  std::cout<<strain_energy_diff_Gibbs<<std::endl;
  std::cout<<strain_energy_Gibbs_new<<std::endl;
  std::cout<<strain_energy_diff_Helm<<std::endl;
  std::cout<<strain_energy_Helm_new<<std::endl;







  

  return EXIT_SUCCESS;
}

// EOF
