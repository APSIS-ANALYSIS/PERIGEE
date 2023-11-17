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
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "QuadPts_debug.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_quad4.hpp"
#include "QuadPts_vis_quad9.hpp"
#include "QuadPts_Gauss_Tet.hpp"

int main(int argc, char *argv[])
{
  QuadPts_vis_quad4 test4{};
  QuadPts_vis_quad9 test9{};

  test4.print_info();
  test9.print_info();

  QuadPts_Gauss_Tet a(5);

  a.print_info();

  QuadPts_vis_tet4 b{};
  b.print_info();
  return EXIT_SUCCESS;
}

// EOF
