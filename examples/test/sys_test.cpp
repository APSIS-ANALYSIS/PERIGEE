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

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);

  double xx = 0;
  double yy = 0;
  double zz = 0;

  double array_x[3];
  double array_y[3];
  double array_z[3];
  
  for(int ii=0; ii<3; ++ii)
  {
    array_x[ii] = MATH_T::gen_double_rand(-1.0, 1.0);
  }

  for(int ii=0; ii<3; ++ii)
  {
    array_y[ii] = MATH_T::gen_double_rand(-1.0, 1.0);
  }

  for(int ii=0; ii<3; ++ii)
  {
    array_z[ii] = MATH_T::gen_double_rand(-1.0, 1.0);
  }
    
// void 类型
  double r1 = 0;  
  FE_T::get_tet_sphere_info(array_x[0], array_x[1], array_x[2], array_x[3],
                            array_y[0], array_y[1], array_y[2], array_y[3],
                            array_z[0], array_z[1], array_z[2], array_z[3],
                            xx, yy, zz, r1);
  std::cout << "void: r1 = " << r1 << '\n';

// double 类型
  double r2 = 0;
  r2 = FE_T::get_tet_sphere_radius(array_x[0], array_x[1], array_x[2], array_x[3],
                                   array_y[0], array_y[1], array_y[2], array_y[3],
                                   array_z[0], array_z[1], array_z[2], array_z[3]);
  std::cout << "double: r2 = " << r2 << '\n';

// 误差分析
  std::cout << std::fixed << std::setprecision(8) << "r_error = " << r1 - r2 << '\n';

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
