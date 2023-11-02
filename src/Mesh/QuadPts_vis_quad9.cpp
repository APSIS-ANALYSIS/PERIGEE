#include "QuadPts_vis_quad9.hpp"

QuadPts_vis_quad9::QuadPts_vis_quad9()
{
  qp[0 ] = 0.0;  qp[1 ] = 0.0;  
  qp[2 ] = 1.0;  qp[3 ] = 0.0;
  qp[4 ] = 1.0;  qp[5 ] = 1.0;
  qp[6 ] = 0.0;  qp[7 ] = 1.0;
  qp[8 ] = 0.5;  qp[9 ] = 0.0;  
  qp[10] = 1.0;  qp[11] = 0.5;
  qp[12] = 0.5;  qp[13] = 1.0;
  qp[14] = 0.0;  qp[15] = 0.5;
  qp[16] = 0.5;  qp[17] = 0.5;

  qw[0] = 0.5;
  qw[1] = 0.5;
  qw[2] = 0.5;
  qw[3] = 0.5;
  qw[4] = 0.5;
  qw[5] = 0.5;
  qw[6] = 0.5;
  qw[7] = 0.5;
  qw[8] = 0.5;
}

QuadPts_vis_quad9::~QuadPts_vis_quad9()
{}

void QuadPts_vis_quad9::print_info() const
{
   SYS_T::commPrint("\n===== Visualization Points for Quad4 ===== \n");
  for(int ii=0; ii<9; ++ii)
    SYS_T::commPrint("%e, %e, %e \n", qw[ii], qp[2*ii], qp[2*ii+1]);
  SYS_T::commPrint("========================================= \n");
}

// EOF
