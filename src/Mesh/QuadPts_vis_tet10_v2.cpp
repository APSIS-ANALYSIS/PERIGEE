#include "QuadPts_vis_tet10_v2.hpp"

QuadPts_vis_tet10_v2::QuadPts_vis_tet10_v2()
{
  qp[0]  = 0.0; qp[1]  = 0.0; qp[2]  = 0.0; qp[3]  = 1.0;
  qp[4]  = 1.0; qp[5]  = 0.0; qp[6]  = 0.0; qp[7]  = 0.0;
  qp[8]  = 0.0; qp[9]  = 1.0; qp[10] = 0.0; qp[11] = 0.0;
  qp[12] = 0.0; qp[13] = 0.0; qp[14] = 1.0; qp[15] = 0.0;
  
  qp[16] = 0.5; qp[17] = 0.0; qp[18] = 0.0; qp[19] = 0.5; // 4
  qp[20] = 0.5; qp[21] = 0.5; qp[22] = 0.0; qp[23] = 0.0; // 5
  qp[24] = 0.0; qp[25] = 0.5; qp[26] = 0.0; qp[27] = 0.5; // 6
  qp[28] = 0.0; qp[29] = 0.0; qp[30] = 0.5; qp[31] = 0.5; // 7
  qp[32] = 0.5; qp[33] = 0.0; qp[34] = 0.5; qp[35] = 0.0; // 8
  qp[36] = 0.0; qp[37] = 0.5; qp[38] = 0.5; qp[39] = 0.0; // 9

  const double weight = 0.1 / 6.0;
  for(int ii=0; ii<10; ++ii) qw[ii] = weight;
}

void QuadPts_vis_tet10_v2::print_info() const
{
  SYS_T::commPrint("\n===== Visualization Points for Tet10 v2 ===== \n");
  for(int ii=0; ii<10; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e, %e, %e, %e \n",
        qw[ii], qp[4*ii], qp[4*ii+1], qp[4*ii+2], qp[4*ii+3]);
  SYS_T::commPrint("========================================= \n");
}

// EOF
