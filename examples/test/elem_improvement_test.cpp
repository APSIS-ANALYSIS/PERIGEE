#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "FEAElement_Line2_3D_der0.hpp"
#include "FEAElement_Line3_3D_der0.hpp"
#include "FE_Tools.hpp"
#include "QuadPts_debug.hpp"

int main( int argc, char * argv[] )
{ 
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  double r = 0.237;
  double s = 0.15;
  IQuadPts * quad = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0} );
  /*
  //FEAElement * tri6_3d = new FEAElement_Triangle3_3D_der0( 1 );
  //FEAElement * tri6_3d = new FEAElement_Triangle3_membrane( 1 );
  //FEAElement * tri6_3d = new FEAElement_Triangle6_3D_der0( 1 );
  FEAElement * tri6_3d = new FEAElement_Triangle6_membrane( 1 );
  
  std::vector<double> tri6_x(6), tri6_y(6), tri6_z(6);
  srand(12);
  for(int ii=0; ii<6; ++ii)
  {
    tri6_x[ii] = rand();
    tri6_y[ii] = rand();
    tri6_z[ii] = rand();
    std::cout<<tri6_x[ii]<<"\t";
    std::cout<<tri6_y[ii]<<"\t";
    std::cout<<tri6_z[ii]<<"\t";
    std::cout<<std::endl;
  }
  //std::vector<double> tri6_x {0.5, 1.0, 0.0, 0.5, 0.5, 0.2};
  //std::vector<double> tri6_y {0.0, 0.3, 1.0, 0.2, 0.5, 0.5};
  //std::vector<double> tri6_z {0.0, 0.0, -0.1, 0.3, 0.0, 0.0};
  
  tri6_3d -> buildBasis( quad, &tri6_x[0], &tri6_y[0], &tri6_z[0] );
  std::cout << "==build basis==" << std::endl;

  std::vector<double> R = tri6_3d -> get_R(0);
  std::cout << "R:" << std::endl;
  for(int ii=0; ii<VEC_T::get_size(R); ++ii) std::cout << std::fixed << std::setprecision(16) << R[ii] << std::endl;

  double area = 0.0;
  const Vector_3 sur_pt( tri6_x[0], tri6_y[0], tri6_z[0] );
  const Vector_3 int_pt( 1.0, 1.0, 1.0 );
  Vector_3 out_n = tri6_3d -> get_normal_out(0, sur_pt, int_pt, area);
  
  std::cout << "out normal:" << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(0) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(1) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(2) << std::endl;

  Tensor2_3D Q = tri6_3d -> get_rotationMatrix(0);
  std::cout << "Q:" << std::endl;
  Q.print();
  */

  // test for line elements
  FEAElement * line_elem = new FEAElement_Line3_3D_der0( 1 );
  std::vector<double> ctrl_x(3), ctrl_y(3), ctrl_z(3);
  srand(1);
  for(int ii=0; ii<3; ++ii)
  {
    ctrl_x[ii] = rand();
    ctrl_y[ii] = rand();
    ctrl_z[ii] = rand();
    std::cout<<ctrl_x[ii]<<"\t";
    std::cout<<ctrl_y[ii]<<"\t";
    std::cout<<ctrl_z[ii]<<"\t";
    std::cout<<std::endl;
  }

  line_elem -> buildBasis( quad, &ctrl_x[0], &ctrl_y[0], &ctrl_z[0] );
  std::cout << "==build basis==" << std::endl;
  
  double * const R = new double [2];
  line_elem -> get_R(0, R);
  std::cout << "R:" << std::endl;
  for(int ii=0; ii<2; ++ii) std::cout << std::fixed << std::setprecision(16) << R[ii] << std::endl;

  double lenth = 0.0;
  std::vector<Vector_3> ctrl_pt { Vector_3(ctrl_x[0], ctrl_y[0], ctrl_z[0]),
	                          Vector_3(ctrl_x[1], ctrl_y[1], ctrl_z[1]),
				  Vector_3(ctrl_x[2], ctrl_y[2], ctrl_z[2]) };

  const Vector_3 int_pt( 0.0, 0.0, 0.0 );
  Vector_3 out_n = line_elem -> get_normal_out( 0, ctrl_pt, int_pt, lenth );

  std::cout << "out normal:" << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(0) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(1) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(2) << std::endl;

  delete line_elem; delete[] R; delete quad;
  PetscFinalize();
  return EXIT_SUCCESS;
}
