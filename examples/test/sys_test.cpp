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
#include "MaterialModel_vol_M94.hpp"
#include "MaterialModel_GOH06_Incompressible_Mixed.hpp"
#include "MaterialModel_GOH14_ST91_Mixed.hpp"
#include "MaterialModel_ich_GOH06.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "MaterialModel_ich_StVenant_Kirchhoff.hpp"
#include <memory>

int main(int argc, char *argv[])
{
  const double rho0 = 1.52333;
  const double elastic_E = 1.2533e1; 
  const double elastic_nu = 0.2832;

  const double f1_the=52.28, f1_phi=-10.98, f2_the=-39.98, f2_phi=-69.98;
  const double fk1=1.6, fk2=0.533, fkd=0.1225; 

  IMaterialModel * oldmodel = new MaterialModel_GOH14_ST91_Mixed(rho0, elastic_E, elastic_nu, f1_the, f1_phi, f2_the, f2_phi, fk1, fk2, fkd);

  std::unique_ptr<IMaterialModel_ich> imodel = SYS_T::make_unique<MaterialModel_ich_GOH14>(oldmodel->get_elastic_mu(), f1_the, f1_phi, f2_the, f2_phi, fk1, fk2, fkd);

  std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_ST91>(rho0, oldmodel->get_elastic_kappa());

  MaterialModel_Mixed_Elasticity * matmodel = new MaterialModel_Mixed_Elasticity(std::move(vmodel), std::move(imodel));
  
  matmodel->print_info();

  std::cout<<std::endl;
  oldmodel->print_info();

  std::cout<<matmodel->get_model_name()<<'\n';

  //oldmodel->print_info();
  
  matmodel->get_fibre_dir(0).print();
  matmodel->get_fibre_dir(1).print();

  Tensor2_3D F; F.gen_rand();

  if(F.det() < 0.0) F*=-1.0;

  std::cout<<"F: \t"; 
  F.print_in_row();

  const double val = std::pow(F.det(), -1.0/3.0);

  std::cout<<"detF = \t"<<val<<"\n\n";

  //F *= val; 
  
  Tensor2_3D P_old, S_old;
  oldmodel->get_PK(F, P_old, S_old);
  // S
  //S_old -= S_new;
  auto S_new = matmodel->get_PK_2nd(F);
  S_old -= S_new.full();
  std::cout<<"PK_2d diff: \t"<<std::sqrt(S_old.MatContraction(S_old))<<'\n';
  std::cout<<"PK_2d valu: \t"<<std::sqrt(S_new.MatContraction(S_new))<<"\n\n";

  //S_old.print_in_row();

  // Cauchy
  Tensor2_3D sigma_old = oldmodel->get_Cauchy_stress(F);
  
  auto sigma_new = matmodel->get_Cauchy_stress(F);
  sigma_old -= sigma_new.full();
  std::cout<<"cauchy diff: \t"<<std::sqrt(sigma_old.MatContraction(sigma_old))<<'\n';
  std::cout<<"cauchy valu: \t"<<std::sqrt(sigma_new.MatContraction(sigma_new))<<"\n\n";

  //P_old -= P_new
  auto P_new = matmodel->get_PK_1st(F);
  P_old -= P_new;
  std::cout<<"PK_1t diff: \t"<<std::sqrt(P_old.MatContraction(P_old))<<'\n';
  std::cout<<"PK_1t valu: \t"<<std::sqrt(P_new.MatContraction(P_new))<<'\n';

  //CC_old -= CC_new
  std::cout<<"\nget_PK_Stiffness:"<<std::endl; 
  Tensor4_3D CC_old;
  oldmodel->get_PK_Stiffness(F, P_old, S_old, CC_old);
  auto CC_new = matmodel->get_PK_Stiffness(F, P_new);
  CC_old -= CC_new.full();
  std::cout<<"CC diff:\t"<<std::sqrt(CC_old.Ten4Contraction(CC_old))<<'\n';
  std::cout<<"CC new: \t"<<std::sqrt(CC_new.full().Ten4Contraction(CC_new.full()))<<"\n\n";
  P_old -= P_new;
  std::cout<<"PK_1s diff: \t"<<std::sqrt(P_old.MatContraction(P_old))<<'\n';
  std::cout<<"PK_1s valu: \t"<<std::sqrt(P_new.MatContraction(P_new))<<'\n';

  std::cout<<"\nget_PK_FFStiffness:"<<std::endl; 
  oldmodel->get_PK_FFStiffness(F, P_old, S_old, CC_old);
  auto CC_n = matmodel->get_PK_FFStiffness(F, P_new);
  CC_old -= CC_n;
  //CC_old.print_in_mat();
  //CC_n.print_in_mat();
  //std::cout<<std::endl; 
  std::cout<<"CC diff: \t"<<std::sqrt(CC_old.Ten4Contraction(CC_old))<<'\n';
  std::cout<<"CC new: \t"<<std::sqrt(CC_n.Ten4Contraction(CC_n))<<"\n\n";
  P_old -= P_new;
  std::cout<<"PK_1st diff: \t"<<std::sqrt(P_old.MatContraction(P_old))<<'\n';
  std::cout<<"PK_1st valu: \t"<<std::sqrt(P_new.MatContraction(P_new))<<'\n';
  std::cout<<std::endl;
  S_old -= S_new.full();
  std::cout<<"PK_2d diff: \t"<<std::sqrt(S_old.MatContraction(S_old))<<'\n';
  std::cout<<"PK_2d valu: \t"<<std::sqrt(S_new.MatContraction(S_new))<<'\n';

  double p = MATH_T::gen_double_rand(-1.0,1.0);
  // rho
  std::cout<<"\nget_rho: "; 
  auto rho_old = oldmodel->get_rho(p);
  auto rho_new = matmodel->get_rho(p);
  rho_old -= rho_new;
  std::cout<<'\t'<<rho_old<<'\t'<<rho_new<<std::endl;

  // beta
  std::cout<<"get_beta: "; 
  auto beta_old = oldmodel->get_beta(p);
  auto beta_new = matmodel->get_beta(p);
  beta_old -= beta_new;
  std::cout<<'\t'<<beta_old<<'\t'<<beta_new<<std::endl;

  // drho_dp
  std::cout<<"get_drho_dp: "; 
  auto drho_dp_old = oldmodel->get_drho_dp(p);
  auto drho_dp_new = matmodel->get_drho_dp(p);
  drho_dp_old -= drho_dp_new;
  std::cout<<'\t'<<drho_dp_old<<'\t'<<drho_dp_new<<std::endl;

  // dbeta_dp
  std::cout<<"get_dbeta_dp: "; 
  auto dbeta_dp_old = oldmodel->get_dbeta_dp(p);
  auto dbeta_dp_new = matmodel->get_dbeta_dp(p);
  dbeta_dp_old -= dbeta_dp_new;
  std::cout<<'\t'<<dbeta_dp_old<<'\t'<<dbeta_dp_new<<std::endl;

  // strain_energy
  std::cout<<"get_strain_energy: "; 
  auto strain_energy_old = oldmodel->get_strain_energy(F);
  auto strain_energy_new = matmodel->get_ich_energy(F);
  strain_energy_old -= strain_energy_new;
  std::cout<<'\t'<<strain_energy_old<<'\t'<<strain_energy_new<<std::endl;

  return EXIT_SUCCESS;
}

// EOF
