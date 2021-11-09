#include "FEAElement_Triangle6.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "QuadPts_debug.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  /*
  double r = 0.237;
  double s = 0.15;
  IQuadPts * quadv = new QuadPts_debug( 4, 1, {r, s, 0.0, 1.0-r-s-0.0}, {1.0} );
  IQuadPts * quads = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0});

  FEAElement * elementv = new FEAElement_Tet10_v2( 1 );
  FEAElement * elements = new FEAElement_Triangle6( 1 );

  std::vector<double> ept_x {0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0};
  std::vector<double> ept_y {0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5};
  std::vector<double> ept_z {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5};

  elementv -> buildBasis( quadv, &ept_x[0], &ept_y[0], &ept_z[0] );

  std::vector<double> spt_x { ept_x[0], ept_x[1], ept_x[2], ept_x[4], ept_x[5], ept_x[6] };
  std::vector<double> spt_y { ept_y[0], ept_y[1], ept_y[2], ept_y[4], ept_y[5], ept_y[6] };

  elements -> buildBasis( quads, &spt_x[0], &spt_y[0] );

  VEC_T::print( elementv->get_d2R_dxy(0) );
  VEC_T::print( elements->get_d2R_dxy(0) );
  */
  
  // Test tet4
  double r = 0.228;
  double s = 0.173;
  double t = 0.119;
  IQuadPts * quad4 = new QuadPts_debug( 4, 1, {r, s, t, 1.0-r-s-t}, {1.0} );
  
  FEAElement * tet4 = new FEAElement_Tet4( 1 );
  
  std::vector<double> tet4_x { 0.0, 2.0, -0.90, 0.1 }; 
  std::vector<double> tet4_y { 0.0, 0.0, 0.3, 0.2 }; 
  std::vector<double> tet4_z { 0.0, 0.0, 2.0, 1.375 }; 

  tet4 -> buildBasis( quad4, &tet4_x[0], &tet4_y[0], &tet4_z[0] );

  auto jac = tet4 -> get_Jacobian(0);
  auto inv_jac = tet4 -> get_invJacobian(0);

  double * old_jac = new double [9];
  double * old_inv_jac = new double [9];
  
  tet4 -> get_Jacobian(0, old_jac);
  tet4 -> get_invJacobian(0, old_inv_jac);


  for(int ii=0; ii<9; ++ii)
  {
    std::cout<<jac[ii] - old_jac[ii]<<'\n';
    std::cout<<inv_jac[ii] - old_inv_jac[ii]<<'\n';
  }

  

  delete tet4; delete quad4;
  // end test tet4

  //delete quadv; delete quads; delete elementv; delete elements;
  PetscFinalize();
  return EXIT_SUCCESS;
}
