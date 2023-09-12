#include "FEAElement_Hex8.hpp"
#include "FEAElement_Quad4.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "Vector_3.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // Test quadrature 1D
  const int poly_degree = 5;
  const int num_pts = 4;
  const double lower_bound = 0;
  const double upper_bound = 1.0;
  double integral = 0.0;
  double integral_qua = 0.0;

  QuadPts_Gauss quad_gauss(num_pts, lower_bound, upper_bound);

  double * coefficient = new double[6] {1.0, 6.9, 3.14, 0.1, 0.05, 0.34};
  int count = 0;
  for(int ii = 0; ii<=poly_degree; ++ii)
  {
    double temp = coefficient[ii];
    for(int qua = 0; qua<num_pts; ++qua)
    {
      integral_qua += temp * std::pow(quad_gauss.get_qp(qua), ii) * quad_gauss.get_qw(qua);  
    }
    integral += (std::pow(upper_bound, ii+1) - std::pow(lower_bound, ii+1)) * temp/(ii+1);
  }
  std::cout <<"Testing Quadrature 1D:" << std::endl;
  std::cout << std::setprecision(16) << integral_qua << std::endl;
  std::cout << std::setprecision(16) << integral << std::endl;
  std::cout << std::setprecision(16) << integral_qua - integral << std::endl;
  
  // Test quadrature 2D
  // Quad
  const int poly_degree_x = 2;
  const int poly_degree_y = 2;
  const int num_pts_x = 2;
  const int num_pts_y = 2;
  const double lower_bound_x = 0.0;
  const double upper_bound_x = 1.0;
  const double lower_bound_y = -1.0;
  const double upper_bound_y = 0.0;
  double integral_quad = 0.0;
  double integral_qua_quad = 0.0;

  QuadPts_Gauss_Quad quad_gauss_quad(num_pts_x, num_pts_y,
    lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y);
  
  count = 0;
  for(int ii = 0; ii<=poly_degree_x; ++ii)
  {
    for(int jj = 0; jj<=poly_degree_y; ++jj)
    {
      double temp = coefficient[count];
      ++count;
      for(int qua = 0; qua<num_pts_x*num_pts_y; ++qua)
      {
        integral_qua_quad += temp * std::pow(quad_gauss_quad.get_qp(qua, 0), ii)
                                  * std::pow(quad_gauss_quad.get_qp(qua, 1), jj)
                                  * quad_gauss_quad.get_qw(qua);
      }
      integral_quad += (std::pow(upper_bound_x, ii+1) - std::pow(lower_bound_x, ii+1)) / (ii+1)
                     * (std::pow(upper_bound_y, jj+1) - std::pow(lower_bound_y, jj+1)) / (jj+1)
                     * temp;
    }
  }
  std::cout <<"Testing Quadrature Quad:" << std::endl;
  std::cout << std::setprecision(16) << integral_qua_quad << std::endl;
  std::cout << std::setprecision(16) << integral_quad << std::endl;
  std::cout << std::setprecision(16) << integral_qua_quad - integral_quad << std::endl;

  // Test quadrature 3D
  // Hex
  const int poly_degree_xx = 2;
  const int poly_degree_yy = 2;
  const int poly_degree_zz = 2;
  const int num_pts_xx = 2;
  const int num_pts_yy = 2;
  const int num_pts_zz = 2;
  const double lower_bound_xx = 0.0;
  const double upper_bound_xx = 1.0;
  const double lower_bound_yy = 0.0;
  const double upper_bound_yy = 1.0;
  const double lower_bound_zz = 0.0;
  const double upper_bound_zz = 1.0;
  double integral_hex = 0.0;
  double integral_qua_hex = 0.0;

  QuadPts_Gauss_Hex quad_gauss_hex(num_pts_xx, num_pts_yy, num_pts_zz,
    lower_bound_xx, upper_bound_xx, lower_bound_yy, upper_bound_yy,
    lower_bound_zz, upper_bound_zz);
  count = 0;
  for(int ii = 0; ii<=poly_degree_xx; ++ii)
  {
    for(int jj = 0; jj<=poly_degree_yy; ++jj)
    {
      for(int kk = 0; kk<=poly_degree_zz; ++kk)
      {
        double temp = coefficient[count];
        ++count;
        for(int qua = 0; qua<num_pts_xx*num_pts_yy*num_pts_zz; ++qua)
        {
          integral_qua_hex += temp * std::pow(quad_gauss_hex.get_qp(qua, 0), ii)
                                   * std::pow(quad_gauss_hex.get_qp(qua, 1), jj)
                                   * std::pow(quad_gauss_hex.get_qp(qua, 2), kk)
                                   * quad_gauss_hex.get_qw(qua);
        }
        integral_hex += (std::pow(upper_bound_xx, ii+1) - std::pow(lower_bound_xx, ii+1)) / (ii+1)
                      * (std::pow(upper_bound_yy, jj+1) - std::pow(lower_bound_yy, jj+1)) / (jj+1)
                      * (std::pow(upper_bound_zz, kk+1) - std::pow(lower_bound_zz, kk+1)) / (kk+1)
                      * temp;
      }
    }
  }
  std::cout <<"Testing Quadrature 3D Hex:" << std::endl;
  std::cout << std::setprecision(16) << integral_qua_hex << std::endl;
  std::cout << std::setprecision(16) << integral_hex << std::endl;
  std::cout << std::setprecision(16) << integral_qua_hex - integral_hex << std::endl;

  // Test Hex element
  Vector_3 vec1,vec2,vec3;
  vec1.gen_rand();
  vec2.gen_rand();
  vec3.gen_rand();

  Vector_3 vec_cross = cross_product(vec1, vec2);
  if(vec_cross.dot_product( vec3 )<0.0)
  {
    Vector_3 temp = vec1;
    vec1 = vec2;
    vec2 = temp;
  }

  const double volume = vec3.dot_product(cross_product(vec1, vec2));

  double * ctrl_x_hex = new double[8]{0.0, vec1.x(), vec1.x()+vec2.x(), vec2.x(), 
    vec3.x(), vec3.x()+vec1.x(), vec3.x()+vec1.x()+vec2.x(), vec3.x()+vec2.x()};
  double * ctrl_y_hex = new double[8]{0.0, vec1.y(), vec1.y()+vec2.y(), vec2.y(), 
    vec3.y(), vec3.y()+vec1.y(), vec3.y()+vec1.y()+vec2.y(), vec3.y()+vec2.y()};
  double * ctrl_z_hex = new double[8]{0.0, vec1.z(), vec1.z()+vec2.z(), vec2.z(), 
    vec3.z(), vec3.z()+vec1.z(), vec3.z()+vec1.z()+vec2.z(), vec3.z()+vec2.z()};
  
  FEAElement_Hex8 hex1(1);
  IQuadPts * quad1 = new QuadPts_Gauss_Hex(1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  hex1.buildBasis(quad1, ctrl_x_hex, ctrl_y_hex, ctrl_z_hex);
  double volume1 = 0.0;
  volume1 += quad1->get_qw(0) * hex1.get_detJac(0);

  FEAElement_Hex8 hex2(8);
  IQuadPts * quad2 = new QuadPts_Gauss_Hex(2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  hex2.buildBasis(quad2, ctrl_x_hex, ctrl_y_hex, ctrl_z_hex);
  double volume2 = 0.0;
  double * sum2 = new double[8]{};
  for(int qua = 0; qua<8; ++qua)
  {
    volume2 += quad2->get_qw(qua) * hex2.get_detJac(qua);
    auto base = hex2.get_R(qua);
    for(int nbas = 0; nbas<8; ++nbas) sum2[qua] += base[nbas]; 
  }

  FEAElement_Hex8 hex3(27);
  IQuadPts * quad3 = new QuadPts_Gauss_Hex(3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  hex3.buildBasis(quad3, ctrl_x_hex, ctrl_y_hex, ctrl_z_hex);
  double volume3 = 0.0;
  double * sum3 = new double[27]{};
  for(int qua = 0; qua<27; ++qua)
  {
    volume3 += quad3->get_qw(qua) * hex3.get_detJac(qua);
    auto base = hex3.get_R(qua);
    for(int nbas = 0; nbas<8; ++nbas) sum3[qua] += base[nbas];
  }

  double * ctrl_x_quad_3D = new double[4]{0.0, vec1.x(), vec1.x()+vec2.x(), vec2.x()};
  double * ctrl_y_quad_3D = new double[4]{0.0, vec1.y(), vec1.y()+vec2.y(), vec2.y()};
  double * ctrl_z_quad_3D = new double[4]{0.0, vec1.z(), vec1.z()+vec2.z(), vec2.z()};
  FEAElement_Quad4_3D_der0 quad_3D(4);
  IQuadPts * quad_q = new QuadPts_Gauss_Quad(2, 0.0, 1.0, 0.0, 1.0);
  quad_3D.buildBasis(quad_q, ctrl_x_quad_3D, ctrl_y_quad_3D, ctrl_z_quad_3D);
  double * sum_quad_3D = new double[4]{};
  for(int qua = 0; qua<4; ++qua)
  {
    auto base = quad_3D.get_R(qua);
    for(int nbas = 0; nbas<4; ++nbas) sum_quad_3D[qua] += base[nbas];
  }

  std::cout << "Testing hex element:" << std::endl;
  hex2.print_info();
  // test buildBasis
  // test R by the sum of basis
  for(int ii = 0; ii<8; ii++)  std::cout << std::setprecision(16) << sum2[ii] << " ";
  std::cout << std::endl;
  for(int ii = 0; ii<27; ii++)  std::cout << std::setprecision(16) << sum3[ii] << " ";
  std::cout << std::endl;
  // test the first order derivative of basis by evaluate the volume of a hex element
  std::cout << std::setprecision(16) << volume - volume1 << std::endl;
  std::cout << std::setprecision(16) << volume - volume2 << std::endl;
  std::cout << std::setprecision(16) << volume - volume3 << std::endl;

  std::cout << "Testing quad 3D element:" << std::endl;
  // test buildBasis
  // test R by the sum of basis
  for(int ii = 0; ii<4; ii++)  std::cout << std::setprecision(16) << sum_quad_3D[ii] << " ";
  std::cout << std::endl;
  // test get_2d_normal_out
  Vector_3 norm_vec = cross_product(vec1,vec2);
  norm_vec.normalize();
  double area = 1.0;
  Vector_3 test_norm_vec = quad_3D.get_2d_normal_out(1, area);
  norm_vec.print();
  test_norm_vec.print();

  PetscFinalize();

  delete quad1; delete quad2; delete quad3;
  delete quad_q;
  delete[] coefficient;
  delete[] ctrl_x_hex;
  delete[] ctrl_y_hex;
  delete[] ctrl_z_hex;
  delete[] ctrl_x_quad_3D;
  delete[] ctrl_y_quad_3D;
  delete[] ctrl_z_quad_3D;
  delete[] sum2; delete[] sum3;
  delete[] sum_quad_3D;
  return EXIT_SUCCESS;
}
