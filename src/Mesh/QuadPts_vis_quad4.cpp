#include "QuadPts_vis_quad4.hpp"

QuadPts_vis_quad4::QuadPts_vis_quad4()
: num_pts( 2 )
{
  qp.resize( 2 * num_pts * num_pts );
  qw.resize( num_pts * num_pts );

  qp[0] = 0.0;  qp[1] = 0.0;  
  qp[2] = 1.0;  qp[3] = 0.0;
  qp[4] = 1.0;  qp[5] = 1.0;
  qp[6] = 0.0;  qp[7] = 1.0;

  qw[0] = 0.5;
  qw[1] = 0.5;
  qw[2] = 0.5;
  qw[3] = 0.5;
}


QuadPts_vis_quad4::~QuadPts_vis_quad4()
{
}


void QuadPts_vis_quad4::print_info() const
{
   SYS_T::commPrint("\n===== Visualization Points for Quad4 ===== \n");
  for(int ii=0; ii<4; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e, %e \n",
        qw[ii], qp[2*ii], qp[2*ii+1]);
  SYS_T::commPrint("========================================= \n");
}

// EOF
