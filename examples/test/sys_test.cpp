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

void test(const double * const &vv, const unsigned int len)
{
  for(unsigned int ii=0; ii<len; ++ii)
    std::cout<<vv[ii]<<'\t';
  std::cout<<"End of vector \n";
}

int main(int argc, char *argv[])
{
  FEAElement * aa = new FEAElement_Tet4( 4 );

  double vx[4] = {1.0, 2.0, 3.0, 4.0};
  double vy[4] = {1.1, 2.1, 3.1, 4.1};
  double vz[4] = {1.2, 2.2, 3.2, 4.3};

  std::vector<double> fx, fy, fz;

  const int id = 2;

  aa -> get_face_ctrlPts( id, &vx[0], &vy[0], &vz[0], fx, fy, fz );

  auto out = aa -> get_face_ctrlPts( id, &vx[0], &vy[0], &vz[0] );

  if( VEC_T::is_equal(fx, out[0], 1.0e-18)) std::cout<<"YES\n";
  else std::cout<<"NO\n"; 

  if( VEC_T::is_equal(fy, out[1], 1.0e-18)) std::cout<<"YES\n";
  else std::cout<<"NO\n"; 
  
  if( VEC_T::is_equal(fz, out[2], 1.0e-18)) std::cout<<"YES\n";
  else std::cout<<"NO\n"; 
 
  test(&out[0][0], out[0].size());
  test(&out[1][0], out[1].size());
  test(&out[2][0], out[2].size());
   
  return EXIT_SUCCESS;
}

// EOF
