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
#include "MaterialModel_ich_GOH06.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "MaterialModel_ich_StVenant_Kirchhoff.hpp"

int main(int argc, char *argv[])
{
  std::vector<double> all;
  all.clear();
  for(int ii=0; ii<10000; ++ii)
  {
    auto v = STen2::gen_rand(-1, 3);
    VEC_T::insert_end(all, v.to_std_vector());
  }

  MATH_T::print_Histogram(all);

  std::cout<<all.size()<<'\n';

  auto AA = Ten4::gen_rand(-8,11);
  auto A = AA;
  auto B = AA;

  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      for(int kk=0; kk<3; ++kk)
        for(int ll=0; ll<3; ++ll)
          B(27*ii+9*jj+3*kk+ll) = 0.5*(AA(27*ii+9*jj+3*kk+ll) + AA(27*ii+9*jj+3*ll+kk));

  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      for(int kk=0; kk<3; ++kk)
        for(int ll=0; ll<3; ++ll)
          A(27*ii+9*jj+3*kk+ll) = 0.5*(B(27*ii+9*jj+3*kk+ll) + B(27*jj+9*ii+3*kk+ll));

  if(A.is_major_sym()) std::cout<<"Yes\n";
  else std::cout<<"No\n";

  if(A.is_minor_sym()) std::cout<<"Yes\n";
  else std::cout<<"No\n";

  auto LHS = STen2::gen_rand(1.0, 2.5);

  auto sol1 = A.solve(LHS.full());
  auto sol2 = A.solve(LHS);

  sol1.print();

  sol1 -= sol2.full();

  sol1.print();
  
  return EXIT_SUCCESS;
}

// EOF
