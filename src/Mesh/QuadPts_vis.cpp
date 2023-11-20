#include "QuadPts_vis.hpp"

QuadPts_vis::QuadPts_vis( const int &input_num_pts )
{
  num_pts = input_num_pts;

  SYS_T::print_fatal_if( (num_pts<2) , "Error: the number of visualization sampling points is less than 2. \n" );

  const double step = 1.0 / (num_pts - 1);

  const double weight = 1.0 / num_pts;

  qp.push_back(0.0);
  qw.push_back(weight);

  for(int ii=1; ii<num_pts; ++ii)
  {
    qp.push_back(ii * step);
    qw.push_back(weight);
  }
  
  VEC_T::shrink2fit(qp);
  VEC_T::shrink2fit(qw);
}

void QuadPts_vis::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\n===== Visualization Points ===== \n");
  for(int ii=0; ii<num_pts; ++ii)
    PetscPrintf(PETSC_COMM_WORLD, "%e, %e \n", qp[ii], qw[ii]);
  PetscPrintf(PETSC_COMM_WORLD, "================================ \n");
}

// EOF
