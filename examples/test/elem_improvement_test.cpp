#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "QuadPts_debug.hpp"

int main( int argc, char * argv[] )
{  
  double r = 0.237;
  double s = 0.15;
  IQuadPts * quad = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0} );

  //FEAElement * tri6_3d = new FEAElement_Triangle3_3D_der0( 1 );
  //FEAElement * tri6_3d = new FEAElement_Triangle3_membrane( 1 );
  //FEAElement * tri6_3d = new FEAElement_Triangle6_3D_der0( 1 );
  FEAElement * tri6_3d = new FEAElement_Triangle6_membrane( 1 );

  std::vector<double> tri6_x(6), tri6_y(6), tri6_z(6);
  srand(5);
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

  std::vector<double> R = tri6_3d -> get_R(0);
  std::cout.precision(16);
  for(int ii=0; ii<R.size(); ++ii) std::cout << std::fixed << std::setprecision(16) << R[ii] << std::endl;

  double area = 0.0;
  const Vector_3 sur_pt( tri6_x[0], tri6_y[0], tri6_z[0] );
  const Vector_3 int_pt( 1.0, 1.0, 1.0 );
  Vector_3 out_n = tri6_3d -> get_normal_out(0, sur_pt, int_pt, area);
  
  std::cout << std::fixed << std::setprecision(16) << out_n(0) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(1) << std::endl;
  std::cout << std::fixed << std::setprecision(16) << out_n(2) << std::endl;

  Tensor2_3D Q = tri6_3d -> get_rotationMatrix(0);
  Q.print();

  delete tri6_3d; delete quad;
  return EXIT_SUCCESS;
}
