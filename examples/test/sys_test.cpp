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
#include "MaterialModel_GOH06_ST91_Mixed.hpp"
#include "MaterialModel_ich_GOH06.hpp"
#include "MaterialModel_ich_GOH14.hpp"
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
  
  std::cout<<"get_PK_2nd:"<<std::endl; 
  Tensor2_3D P_old, S_old;
  oldmodel->get_PK(F, P_old, S_old);
  // S
  //S_old -= S_new;
  auto S_new = matmodel->get_PK_2nd(F);
  S_old -= S_new.full();
  S_old.print_in_row();
  std::cout<<std::endl;
  S_new.print_in_row();
  std::cout<<std::endl; 

  //P_old -= P_new
  std::cout<<"get_PK_1st:"<<std::endl; 
  auto P_new = matmodel->get_PK_1st(F);
  P_old -= P_new;
  P_old.print_in_row();
  std::cout<<std::endl;
  P_new.print_in_row();
  std::cout<<std::endl; 

  //CC_old -= CC_new
  std::cout<<"get_PK_Stiffness:"<<std::endl; 
  Tensor4_3D CC_old;
  oldmodel->get_PK_Stiffness(F, P_old, S_old, CC_old);
  auto CC_new = matmodel->get_PK_Stiffness(F, P_new);
  CC_old -= CC_new.full();
  CC_old.print_in_mat();
  std::cout<<std::endl;
  CC_new.print_in_mat();
  std::cout<<std::endl; 
  P_old -= P_new;
  P_old.print_in_row();
  std::cout<<std::endl;
  P_new.print_in_row();
  std::cout<<std::endl;

  std::cout<<"get_PK_FFStiffness:"<<std::endl; 
  oldmodel->get_PK_FFStiffness(F, P_old, S_old, CC_old);
  auto CC_n = matmodel->get_PK_FFStiffness(F, P_new);
  CC_old -= CC_n;
  CC_old.print_in_mat();
  std::cout<<std::endl;
  CC_n.print_in_mat();
  std::cout<<std::endl; 
  P_old -= P_new;
  P_old.print_in_row();
  std::cout<<std::endl;
  P_new.print_in_row();
  std::cout<<std::endl;
  S_old -= S_new.full();
  S_old.print_in_row();

  double p = MATH_T::gen_double_rand(-1.0,1.0);
  // rho
  std::cout<<"get_rho:"<<std::endl; 
  auto rho_old = oldmodel->get_rho(p);
  auto rho_new = matmodel->get_rho(p);
  rho_old -= rho_new;
  std::cout<<rho_old<<std::endl;
  std::cout<<rho_new<<std::endl;

  // beta
  std::cout<<"get_beta:"<<std::endl; 
  auto beta_old = oldmodel->get_beta(p);
  auto beta_new = matmodel->get_beta(p);
  beta_old -= beta_new;
  std::cout<<beta_old<<std::endl;
  std::cout<<beta_new<<std::endl;

  // drho_dp
  std::cout<<"get_drho_dp:"<<std::endl; 
  auto drho_dp_old = oldmodel->get_drho_dp(p);
  auto drho_dp_new = matmodel->get_drho_dp(p);
  drho_dp_old -= drho_dp_new;
  std::cout<<drho_dp_old<<std::endl;
  std::cout<<drho_dp_new<<std::endl;

  // dbeta_dp
  std::cout<<"get_dbeta_dp:"<<std::endl; 
  auto dbeta_dp_old = oldmodel->get_dbeta_dp(p);
  auto dbeta_dp_new = matmodel->get_dbeta_dp(p);
  dbeta_dp_old -= dbeta_dp_new;
  std::cout<<dbeta_dp_old<<std::endl;
  std::cout<<dbeta_dp_new<<std::endl;

  // strain_energy
  std::cout<<"get_strain_energy:"<<std::endl; 
  auto strain_energy_old = oldmodel->get_strain_energy(F);
  auto strain_energy_new = matmodel->get_ich_energy(F);
  strain_energy_old -= strain_energy_new;
  std::cout<<strain_energy_old<<std::endl;
  std::cout<<strain_energy_new<<std::endl;

  /*
  // GOH06
  const double f1_the=49.98, f1_phi=49.98, f2_the=-49.98, f2_phi=-49.98;
  const double fk1=996.6, fk2=0.1, fkd=1.0 /3.0; 

  IMaterialModel * oldmodel_GOH06_ST91 = new MaterialModel_GOH06_ST91_Mixed(rho0, elastic_E, elastic_nu, f1_the, f1_phi, f2_the, f2_phi, fk1, fk2, fkd);

  std::unique_ptr<IMaterialModel_ich> imodel_GOH06 = SYS_T::make_unique<MaterialModel_ich_GOH06>(oldmodel->get_elastic_mu(), f1_the, f1_phi, f2_the, f2_phi, fk1, fk2, fkd);

  std::unique_ptr<IMaterialModel_vol> vmodel_ST91 = SYS_T::make_unique<MaterialModel_vol_ST91>(rho0, oldmodel->get_elastic_kappa());

  MaterialModel_Mixed_Elasticity * matmodel_GOH06_ST91 = new MaterialModel_Mixed_Elasticity(std::move(vmodel_ST91), std::move(imodel_GOH06));

  matmodel_GOH06_ST91->print_info();

  std::cout<<matmodel_GOH06_ST91->get_model_name()<<'\n';

  //oldmodel->print_info();

  F.gen_rand();

  if(F.det() < 0.0) F*=-1.0;

  std::cout<<"F:"<<std::endl; 
  F.print();
  std::cout<<std::endl; 

  // F *= val; 

  P_old.gen_id();
  S_old.gen_id();
  oldmodel_GOH06_ST91->get_PK(F, P_old, S_old);
  // S
  //S_old -= S_new;
  S_new = matmodel_GOH06_ST91->get_PK_2nd(F);
  S_old -= S_new.full();
  S_old.print();
  std::cout<<std::endl;
  S_new.print();
  std::cout<<std::endl; 


  //P_old -= P_new
  P_new = matmodel_GOH06_ST91->get_PK_1st(F);
  P_old -= P_new;
  P_old.print();
  std::cout<<std::endl;
  P_new.print();
  std::cout<<std::endl; 

  //CC_old -= CC_new
  oldmodel_GOH06_ST91->get_PK_Stiffness(F, P_old, S_old, CC_old);
  CC_new = matmodel_GOH06_ST91->get_PK_Stiffness(F, P_new);
  CC_old -= CC_new.full();
  CC_old.print();
  std::cout<<std::endl;
  CC_new.print();
  std::cout<<std::endl;

  p = MATH_T::gen_double_rand(-1.0,1.0);
  // rho
  rho_old = oldmodel_GOH06_ST91->get_rho(p);
  rho_new = matmodel_GOH06_ST91->get_rho(p);
  rho_old -= rho_new;
  std::cout<<rho_old<<std::endl;
  std::cout<<rho_new<<std::endl;

  // beta
  beta_old = oldmodel_GOH06_ST91->get_beta(p);
  beta_new = matmodel_GOH06_ST91->get_beta(p);
  beta_old -= beta_new;
  std::cout<<beta_old<<std::endl;
  std::cout<<beta_new<<std::endl;

  // drho_dp
  drho_dp_old = oldmodel_GOH06_ST91->get_drho_dp(p);
  drho_dp_new = matmodel_GOH06_ST91->get_drho_dp(p);
  drho_dp_old -= drho_dp_new;
  std::cout<<drho_dp_old<<std::endl;
  std::cout<<drho_dp_new<<std::endl;

  // dbeta_dp
  dbeta_dp_old = oldmodel_GOH06_ST91->get_dbeta_dp(p);
  dbeta_dp_new = matmodel_GOH06_ST91->get_dbeta_dp(p);
  dbeta_dp_old -= dbeta_dp_new;
  std::cout<<dbeta_dp_old<<std::endl;
  std::cout<<dbeta_dp_new<<std::endl;

  // strain_energy
  strain_energy_old = oldmodel_GOH06_ST91->get_strain_energy(F);
  strain_energy_new = matmodel_GOH06_ST91->get_ich_energy(F);
  strain_energy_old -= strain_energy_new;
  std::cout<<strain_energy_old<<std::endl;
  std::cout<<strain_energy_new<<std::endl;

  // fibre_dir
  const int dir = MATH_T::gen_int_rand(0,1);
  auto a_old = oldmodel_GOH06_ST91->get_fibre_dir(dir);
  auto a_new = matmodel_GOH06_ST91->get_fibre_dir(dir);
  a_old -= a_new;
  a_old.print();
  a_new.print();
  */

  return EXIT_SUCCESS;
}

// EOF
