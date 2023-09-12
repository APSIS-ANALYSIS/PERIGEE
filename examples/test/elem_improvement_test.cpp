#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "QuadPts_debug.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  
  double r = 0.237;
  double s = 0.15;
  IQuadPts * quad = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0} );

  //FEAElement * tri6_3d = new FEAElement_Triangle6_3D_der0( 1 );
  FEAElement * tri6_3d = new FEAElement_Triangle6_membrane( 1 );

  std::vector<double> tri6_x {0.0, 1.0, 0.0, 0.5, 0.5, 0.0};
  std::vector<double> tri6_y {0.0, 0.0, 1.0, 0.0, 0.5, 0.5};
  std::vector<double> tri6_z {0.0, 0.0, 0.0, 0.3, 0.0, 0.0};
  
  tri6_3d -> buildBasis( quad, &tri6_x[0], &tri6_y[0], &tri6_z[0] );

  std::vector<double> R = tri6_3d -> get_R(0);
  for(int ii=0; ii<R.size(); ++ii) std::cout<<R[ii]<<"\t";
  std::cout<<std::endl;

  double area = 0.0;
  const std::vector<double> int_pt { 1.0, 1.0, 1.0 };
  Vector_3 out_n = tri6_3d -> get_normal_out(0, tri6_x[0], tri6_y[0], tri6_z[0],
                  int_pt[0], int_pt[1], int_pt[2], area);

  out_n.print();

  delete tri6_3d; delete quad;
  
  /* 
  double r = 0.237;
  double s = 0.15;
  IQuadPts * quad = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0} );
  
  //FEAElement * tri3_3d = new FEAElement_Triangle3_3D_der0( 1 );
  FEAElement * tri3_3d = new FEAElement_Triangle3_membrane( 1 );
  
  const std::vector<double> tri3_x { 0.0, 2.0, 0.0 }; 
  const std::vector<double> tri3_y { 0.0, 0.0, 0.3 }; 
  const std::vector<double> tri3_z { 0.0, 0.0, 2.0 }; 

  tri3_3d -> buildBasis( quad, &tri3_x[0], &tri3_y[0], &tri3_z[0] );

  std::vector<double> R = tri3_3d -> get_R(0);
  for(int ii=0; ii<R.size(); ++ii) std::cout<<R[ii]<<"\t";
  std::cout<<std::endl;

  double area = 0.0;
  const std::vector<double> sur_pt { 0.0, 2.0, 0.0 };
  Vector_3 out_n = tri3_3d -> get_normal_out(0, tri3_x[0], tri3_y[0], tri3_z[0], 
		  sur_pt[0], sur_pt[1], sur_pt[2], area);

  out_n.print();

  delete tri3_3d; delete quad;
  */

  PetscFinalize();
  return EXIT_SUCCESS;
}
