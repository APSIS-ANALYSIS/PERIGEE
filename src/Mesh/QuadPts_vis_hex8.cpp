#include "QuadPts_vis_hex8.hpp"

QuadPts_vis_hex8::QuadPts_vis_hex8()
{
  qp[0]  = 0.0; qp[1]  = 0.0; qp[2]  = 0.0; 
  qp[3]  = 1.0; qp[4]  = 0.0; qp[5]  = 0.0; 
  qp[6]  = 1.0; qp[7]  = 1.0; qp[8]  = 0.0;
  qp[9]  = 0.0; qp[10] = 1.0; qp[11] = 0.0; 
  qp[12] = 0.0; qp[13] = 0.0; qp[14] = 1.0;
  qp[15] = 1.0; qp[16] = 0.0; qp[17] = 1.0;
  qp[18] = 1.0; qp[19] = 1.0; qp[20] = 1.0;
  qp[21] = 0.0; qp[22] = 1.0; qp[23] = 1.0;

  qw[0] = 0.5; qw[1] = 0.5;
  qw[2] = 0.5; qw[3] = 0.5;
  qw[4] = 0.5; qw[5] = 0.5; 
  qw[6] = 0.5; qw[7] = 0.5;
}

QuadPts_vis_hex8::~QuadPts_vis_hex8()
{}

void QuadPts_vis_hex8::print_info() const
{
  SYS_T::commPrint("\n===== Visualization Points for Hex8 ===== \n");
  for(int ii=0; ii<8; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e, %e, %e \n", 
        qw[ii], qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
  SYS_T::commPrint("========================================= \n");
}

// EOF