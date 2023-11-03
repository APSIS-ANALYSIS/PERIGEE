#include "QuadPts_vis_tri6.hpp"

QuadPts_vis_tri6::QuadPts_vis_tri6()
{
  qp[0] = 0.0;  qp[1] = 0.0;  qp[2] = 1.0;
  qp[3] = 1.0;  qp[4] = 0.0;  qp[5] = 0.0;
  qp[6] = 0.0;  qp[7] = 1.0;  qp[8] = 0.0;
  qp[9] = 0.5;  qp[10] = 0.0; qp[11] = 0.5;
  qp[12] = 0.5; qp[13] = 0.5; qp[14] = 0.0;
  qp[15] = 0.0; qp[16] = 0.5; qp[17] = 0.5;

  qw[0] = 0.5 / 6.0;
  qw[1] = 0.5 / 6.0;
  qw[2] = 0.5 / 6.0;
  qw[3] = 0.5 / 6.0;
  qw[4] = 0.5 / 6.0;
  qw[5] = 0.5 / 6.0;
}

QuadPts_vis_tri6::~QuadPts_vis_tri6()
{
}

void QuadPts_vis_tri6::print_info() const
{
   SYS_T::commPrint("\n===== Visualization Points for Tri6 ===== \n");
  for(int ii=0; ii<6; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e, %e, %e \n",
        qw[ii], qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
  SYS_T::commPrint("========================================= \n");
}

// EOF
