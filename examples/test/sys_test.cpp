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
#include "MaterialModel_mixed_NeoHookean.hpp"
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

  IMaterialModel_mixed * matmodel = new MaterialModel_mixed_NeoHookean(std::move(vmodel), oldmodel->get_elastic_mu());

  //matmodel->print_info();


  //oldmodel->print_info();

  Tensor2_3D F; F.gen_rand();

  if(F.det() < 0.0) F*=-1.0;

  const double val = std::pow(F.det(), -1.0/3.0);

  //std::cout<<val<<'\n';

  F *= val; 

  auto P_new = matmodel->get_PK_1st(F);
  auto S_new = matmodel->get_PK_2nd(F);
  
  Tensor2_3D P_old, S_old;
  oldmodel->get_PK(F, P_old, S_old);

  //P_old -= P_new;
  S_old -= S_new.convert_to_full();

  S_old.print();
  
  std::cout<<std::endl;

  S_new.print();

  return EXIT_SUCCESS;
}

// EOF
