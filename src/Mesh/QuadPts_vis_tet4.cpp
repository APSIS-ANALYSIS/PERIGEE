#include "QuadPts_vis_tet4.hpp"

QuadPts_vis_tet4::QuadPts_vis_tet4()
{
  qp[0]  = 0.0; qp[1]  = 0.0; qp[2]  = 0.0; qp[3]  = 1.0; 
  qp[4]  = 1.0; qp[5]  = 0.0; qp[6]  = 0.0; qp[7]  = 0.0; 
  qp[8]  = 0.0; qp[9]  = 1.0; qp[10] = 0.0; qp[11] = 0.0; 
  qp[12] = 0.0; qp[13] = 0.0; qp[14] = 1.0; qp[15] = 0.0;

  qw[0] = 0.25 / 6.0; qw[1] = 0.25 / 6.0; 
  qw[2] = 0.25 / 6.0; qw[3] = 0.25 / 6.0;
}

void QuadPts_vis_tet4::print_info() const
{
  SYS_T::commPrint("\n===== Visualization Points for Tet4 ===== \n");
  for(int ii=0; ii<4; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e, %e, %e, %e \n", 
        qw[ii], qp[4*ii], qp[4*ii+1], qp[4*ii+2], qp[4*ii+3]);
  SYS_T::commPrint("========================================= \n");
}


// EOF
