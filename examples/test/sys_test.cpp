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
#include "QuadPts_vis_quad4.hpp"

int main(int argc, char *argv[])
{
  QuadPts_vis_quad4 test{};

  test.print_info();

  return EXIT_SUCCESS;
}

// EOF
