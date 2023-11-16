#include "FEAElement_Hex8.hpp"
#include "FEAElement_Quad4.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad9.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_debug.hpp"
#include "Vector_3.hpp"
#include <iostream>
#include <fstream>

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // Test quadrature 1D
  const double left = -1.0;
  const double right = 1.0;
  const int max_num = 3;
  double * coefficient = new double[1024] {};
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  for(int ii = 0; ii<1024; ++ii)
  {
    coefficient[ii] = dis(gen);
  }
  std::cout <<"Testing Quadrature 1D:" << std::endl;
  for(int num_pts = 1; num_pts<max_num*5; ++num_pts)
  {
    const int poly_degree = 2*num_pts-1;
    const double lower_bound = -5.0;
    const double upper_bound = -4.0;
    double integral = 0.0;
    double integral_qua = 0.0;

    QuadPts_Gauss_1D quad_gauss(num_pts, lower_bound, upper_bound);
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
    std::cout <<"Number of quadrature point:" << num_pts << std::endl;
    std::cout << std::setprecision(16) << integral_qua << std::endl;
    std::cout << std::setprecision(16) << integral << std::endl;
    std::cout << std::setprecision(16) << integral_qua - integral << std::endl;
  }
  
  // Test quadrature 2D
  // Quad
  std::cout <<"Testing Quadrature Quad:" << std::endl;
  for(int num_pts_y = 1; num_pts_y<max_num; ++num_pts_y)
  {
    for(int num_pts_x = 1; num_pts_x<max_num; ++num_pts_x)
    {    
      const int poly_degree_x = 2*num_pts_x-1;
      const int poly_degree_y = 2*num_pts_y-1;
      const double lower_bound_x = 0.0;
      const double upper_bound_x = 1.0;
      const double lower_bound_y = -1.0;
      const double upper_bound_y = 0.0;
      double integral_quad = 0.0;
      double integral_qua_quad = 0.0;

      QuadPts_Gauss_Quad quad_gauss_quad(num_pts_x, num_pts_y,
        lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y);
  
      int count = 0;
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
      std::cout <<"Number of quadrature point x:" << num_pts_x << std::endl;
      std::cout <<"Number of quadrature point y:" << num_pts_y << std::endl;
      std::cout << std::setprecision(16) << integral_qua_quad << std::endl;
      std::cout << std::setprecision(16) << integral_quad << std::endl;
      std::cout << std::setprecision(16) << integral_qua_quad - integral_quad << std::endl;
    }
  }

  // Test quadrature 3D
  // Hex
  std::cout <<"Testing Quadrature 3D Hex:" << std::endl;
  for(int num_pts_zz = 1; num_pts_zz<max_num; ++num_pts_zz)
  {
    for(int num_pts_yy = 1; num_pts_yy<max_num; ++num_pts_yy)
    {
      for(int num_pts_xx = 1; num_pts_xx<max_num; ++num_pts_xx)
      {
        const int poly_degree_xx = 2*num_pts_xx-1;
        const int poly_degree_yy = 2*num_pts_yy-1;
        const int poly_degree_zz = 2*num_pts_zz-1;
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
        int count = 0;
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
        std::cout <<"Number of quadrature point x:" << num_pts_xx << std::endl;
        std::cout <<"Number of quadrature point y:" << num_pts_yy << std::endl;
        std::cout <<"Number of quadrature point z:" << num_pts_zz << std::endl;
        std::cout << std::setprecision(16) << integral_qua_hex << std::endl;
        std::cout << std::setprecision(16) << integral_hex << std::endl;
        std::cout << std::setprecision(16) << integral_qua_hex - integral_hex << std::endl;
      }
    }
  }

  // Test Hex element
  Vector_3 vec1,vec2,vec3;
  vec1.gen_rand();
  vec2.gen_rand();
  vec3.gen_rand();
  Vector_3 vec_cross = Vec3::cross_product(vec1, vec2);
  if(vec_cross.dot_product( vec3 )<0.0)
  {
    Vector_3 temp = vec1;
    vec1 = vec2;
    vec2 = temp;
  }

  const double volume = vec3.dot_product(Vec3::cross_product(vec1, vec2));

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

  // test Hex27
  double * ctrl_x_hex_27 = new double[27]{};
  double * ctrl_y_hex_27 = new double[27]{};
  double * ctrl_z_hex_27 = new double[27]{};
  double * R_matlab = new double[27]{};
  double * dR_dx_matlab = new double[27]{};
  double * dR_dy_matlab = new double[27]{};
  double * dR_dz_matlab = new double[27]{};
  double * d2R_dxx_matlab = new double[27]{};
  double * d2R_dyy_matlab = new double[27]{};
  double * d2R_dzz_matlab = new double[27]{};
  double * d2R_dxy_matlab = new double[27]{};
  double * d2R_dxz_matlab = new double[27]{};
  double * d2R_dyz_matlab = new double[27]{};
  double * J_matlab = new double[9]{};
  double * Jinv_matlab = new double[9]{};
  double detJ_matlab = 0.0;
  std::ifstream infile;
  infile.open("../ctrlpts.txt",std::ios::in);
  if (!infile.is_open())
	{
		std::cout << "error" << std::endl;
	}
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> ctrl_x_hex_27[ii];
  }
  for(int jj = 0; jj<27; ++jj)
  {
    infile >> ctrl_y_hex_27[jj];
  }
  for(int kk = 0; kk<27; ++kk)
  {
    infile >> ctrl_z_hex_27[kk];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> R_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> dR_dx_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> dR_dy_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> dR_dz_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dxx_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dyy_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dzz_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dxy_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dxz_matlab[ii];
  }
  for(int ii = 0; ii<27; ++ii)
  {
    infile >> d2R_dyz_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> J_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> Jinv_matlab[ii];
  }
  infile >> detJ_matlab;
  infile.close();
  FEAElement_Hex27 hex_27(1);
  std::vector<double> in_qp{{0.52572, 0.33891, 0.12345}};
  std::vector<double> in_qw{{1}};
  QuadPts_debug * quad_debug = new QuadPts_debug(in_qp, in_qw, 3);
  hex_27.buildBasis(quad_debug, ctrl_x_hex_27, ctrl_y_hex_27, ctrl_z_hex_27);
  double tol = 1e-14;
  bool isSame_R = true;
  bool isSame_dR_dx = true;
  bool isSame_dR_dy = true;
  bool isSame_dR_dz = true;
  bool isSame_d2R_dxx = true;
  bool isSame_d2R_dyy = true;
  bool isSame_d2R_dzz = true;
  bool isSame_d2R_dxy = true;
  bool isSame_d2R_dxz = true;
  bool isSame_d2R_dyz = true;
  bool isSame_J = true;
  bool isSame_Jinv = true;
  bool isSame_detJ = true;

  double * basis = new double[27]{};
  double * dR_dx = new double[27]{};
  double * dR_dy = new double[27]{};
  double * dR_dz = new double[27]{};
  double * d2R_dxx = new double[27]{};
  double * d2R_dyy = new double[27]{};
  double * d2R_dzz = new double[27]{};
  double * d2R_dxy = new double[27]{};
  double * d2R_dxz = new double[27]{};
  double * d2R_dyz = new double[27]{};
  double * J = new double[9]{};
  double * Jinv = new double[9]{};
  double detJ = 0.0;
  hex_27.get_3D_R_dR_d2R(0, basis, dR_dx, dR_dy, dR_dz, d2R_dxx, d2R_dyy, d2R_dzz,
      d2R_dxy, d2R_dxz, d2R_dyz);
  hex_27.get_Jacobian(0, J);
  hex_27.get_invJacobian(0, Jinv);

  for (int ii = 0; ii<27; ++ii)
  {
    if (std::abs(basis[ii]-R_matlab[ii])>tol) isSame_R = false;
    if (std::abs(dR_dx[ii]-dR_dx_matlab[ii])>tol) isSame_dR_dx = false;
    if (std::abs(dR_dy[ii]-dR_dy_matlab[ii])>tol) isSame_dR_dy = false;
    if (std::abs(dR_dz[ii]-dR_dz_matlab[ii])>tol) isSame_dR_dz = false;
    if (std::abs(d2R_dxx[ii]-d2R_dxx_matlab[ii])>tol) isSame_d2R_dxx = false;
    if (std::abs(d2R_dyy[ii]-d2R_dyy_matlab[ii])>tol) isSame_d2R_dyy = false;
    if (std::abs(d2R_dzz[ii]-d2R_dzz_matlab[ii])>tol) isSame_d2R_dzz = false;
    if (std::abs(d2R_dxy[ii]-d2R_dxy_matlab[ii])>tol) isSame_d2R_dxy = false;
    if (std::abs(d2R_dxz[ii]-d2R_dxz_matlab[ii])>tol) isSame_d2R_dxz = false;
    if (std::abs(d2R_dyz[ii]-d2R_dyz_matlab[ii])>tol) isSame_d2R_dyz = false;
  }
  for (int ii = 0; ii<9; ++ii)
  {
    if (std::abs(J[ii]-J_matlab[ii])>tol) isSame_J = false;
    if (std::abs(Jinv[ii]-Jinv_matlab[ii])>tol) isSame_Jinv = false;
  }

  std::cout << "Testing R: " << isSame_R << std::endl; 
  std::cout << "Testing dR_dx: " << isSame_dR_dx << std::endl;
  std::cout << "Testing dR_dy: " << isSame_dR_dy << std::endl; 
  std::cout << "Testing dR_dz: " << isSame_dR_dz << std::endl;
  std::cout << "Testing d2R_dxx: " << isSame_d2R_dxx << std::endl;
  std::cout << "Testing d2R_dyy: " << isSame_d2R_dyy << std::endl; 
  std::cout << "Testing d2R_dzz: " << isSame_d2R_dzz << std::endl;
  std::cout << "Testing d2R_dxy: " << isSame_d2R_dxy << std::endl;
  std::cout << "Testing d2R_dxz: " << isSame_d2R_dxz << std::endl; 
  std::cout << "Testing d2R_dyz: " << isSame_d2R_dyz << std::endl;  
  std::cout << "Testing J: " << isSame_J << std::endl;
  std::cout << "Testing Jinv: " << isSame_Jinv << std::endl;

  // test quad9
  double * ctrl_x_quad_9 = new double[9]{};
  double * ctrl_y_quad_9 = new double[9]{};
  double * R_quad_matlab = new double[9]{};
  double * dR_dx_quad_matlab = new double[9]{};
  double * dR_dy_quad_matlab = new double[9]{};
  double * d2R_dxx_quad_matlab = new double[9]{};
  double * d2R_dyy_quad_matlab = new double[9]{};
  double * d2R_dxy_quad_matlab = new double[9]{};
  double * J_quad_matlab = new double[4]{};
  double * Jinv_quad_matlab = new double[4]{};

  infile.open("/Users/seavegetable/Documents/MATLAB/code_test/ctrlpts_quad.txt",std::ios::in);
  if (!infile.is_open())
	{
		std::cout << "error" << std::endl;
	}
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> ctrl_x_quad_9[ii];
  }
  for(int jj = 0; jj<9; ++jj)
  {
    infile >> ctrl_y_quad_9[jj];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> R_quad_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> dR_dx_quad_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> dR_dy_quad_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> d2R_dxx_quad_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> d2R_dyy_quad_matlab[ii];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> d2R_dxy_quad_matlab[ii];
  }
  for(int ii = 0; ii<4; ++ii)
  {
    infile >> J_quad_matlab[ii];
  }
  for(int ii = 0; ii<4; ++ii)
  {
    infile >> Jinv_quad_matlab[ii];
  }
  infile.close();
  FEAElement_Quad9 quad_9(1);
  std::vector<double> in_qp_quad{{0.2, 0.3}};
  QuadPts_debug * quad_debug_quad = new QuadPts_debug(in_qp_quad, in_qw, 2);
  quad_9.buildBasis(quad_debug_quad, ctrl_x_quad_9, ctrl_y_quad_9);
  isSame_R = true;
  isSame_dR_dx = true;
  isSame_dR_dy = true;
  isSame_d2R_dxx = true;
  isSame_d2R_dyy = true;
  isSame_d2R_dxy = true;
  isSame_J = true;
  isSame_Jinv = true;
  
  double * basis_quad = new double[9]{};
  double * dR_dx_quad = new double[9]{};
  double * dR_dy_quad = new double[9]{};
  double * d2R_dxx_quad = new double[9]{};
  double * d2R_dyy_quad = new double[9]{};
  double * d2R_dxy_quad = new double[9]{};
  double * J_quad = new double[4]{};
  double * Jinv_quad = new double[4]{};
  quad_9.get_2D_R_dR_d2R(0, basis_quad, dR_dx_quad, dR_dy_quad, d2R_dxx_quad, d2R_dyy_quad, d2R_dxy_quad);
  quad_9.get_Jacobian(0, J_quad);
  quad_9.get_invJacobian(0, Jinv_quad);
  tol = 1e-4;
  
  for (int ii = 0; ii<9; ++ii)
  {
    if (std::abs(basis_quad[ii]-R_quad_matlab[ii])>tol) isSame_R = false;
    if (std::abs(dR_dx_quad[ii]-dR_dx_quad_matlab[ii])>tol) isSame_dR_dx = false;
    if (std::abs(dR_dy_quad[ii]-dR_dy_quad_matlab[ii])>tol) isSame_dR_dy = false;
    if (std::abs(d2R_dxx_quad[ii]-d2R_dxx_quad_matlab[ii])>tol) isSame_d2R_dxx = false;
    if (std::abs(d2R_dyy_quad[ii]-d2R_dyy_quad_matlab[ii])>tol) isSame_d2R_dyy = false;
    if (std::abs(d2R_dxy_quad[ii]-d2R_dxy_quad_matlab[ii])>tol) isSame_d2R_dxy = false;
  }
  for (int ii = 0; ii<4; ++ii)
  {
    if (std::abs(J_quad[ii]-J_quad_matlab[ii])>tol) isSame_J = false;
    if (std::abs(Jinv_quad[ii]-Jinv_quad_matlab[ii])>tol) isSame_Jinv = false;
  }
  std::cout << "============================" << std::endl; 
  std::cout << "Testing R: " << isSame_R << std::endl; 
  std::cout << "Testing dR_dx: " << isSame_dR_dx << std::endl;
  std::cout << "Testing dR_dy: " << isSame_dR_dy << std::endl; 
  std::cout << "Testing d2R_dxx: " << isSame_d2R_dxx << std::endl;
  std::cout << "Testing d2R_dyy: " << isSame_d2R_dyy << std::endl; 
  std::cout << "Testing d2R_dxy: " << isSame_d2R_dxy << std::endl;
  std::cout << "Testing J: " << isSame_J << std::endl;
  std::cout << "Testing Jinv: " << isSame_Jinv << std::endl;

  // test quad 9 3D
  double * ctrl_x_quad9_3D = new double[9]{};
  double * ctrl_y_quad9_3D = new double[9]{};
  double * ctrl_z_quad9_3D = new double[9]{};
  double * R_quad_3D_matlab = new double[9]{};
  double * outn_matlab = new double[3]{};
  
  infile.open("/Users/seavegetable/Documents/MATLAB/code_test/ctrlpts_quad_3D.txt",std::ios::in);
  if (!infile.is_open())
	{
		std::cout << "error" << std::endl;
	}
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> ctrl_x_quad9_3D[ii];
  }
  for(int jj = 0; jj<9; ++jj)
  {
    infile >> ctrl_y_quad9_3D[jj];
  }
  for(int jj = 0; jj<9; ++jj)
  {
    infile >> ctrl_z_quad9_3D[jj];
  }
  for(int ii = 0; ii<9; ++ii)
  {
    infile >> R_quad_3D_matlab[ii];
  }
  for(int ii = 0; ii<3; ++ii)
  {
    infile >> outn_matlab[ii];
  }
  infile.close();
  FEAElement_Quad9_3D_der0 quad9_3D(1);
  quad9_3D.buildBasis(quad_debug_quad, ctrl_x_quad_3D, ctrl_y_quad_3D, ctrl_z_quad_3D);
  isSame_R = true;
  bool isSame_outn = true;
  
  double * basis_quad_3D = new double[9]{};
  quad9_3D.get_R(0,basis_quad_3D);
  double area;
  Vector_3 outn = quad9_3D.get_2d_normal_out(0, area);
  double * nout = new double[3]{};
  nout[0] = outn.x();
  nout[1] = outn.y();
  nout[2] = outn.z();
  tol = 1e-14;
  
  for (int ii = 0; ii<9; ++ii)
  {
    if (std::abs(basis_quad_3D[ii]-R_quad_3D_matlab[ii])>tol) isSame_R = false;
    if (std::abs(nout[ii]-outn_matlab[ii])>tol) isSame_dR_dx = false;
  }
  std::cout << "============================" << std::endl; 
  std::cout << "Testing R: " << isSame_R << std::endl; 
  std::cout << "Testing nout: " << isSame_outn << std::endl;

  PetscFinalize();

  delete quad1; delete quad2; delete quad3;
  delete quad_q; delete quad_debug; delete quad_debug_quad;
  delete[] coefficient;
  delete[] ctrl_x_hex;
  delete[] ctrl_y_hex;
  delete[] ctrl_z_hex;
  delete[] ctrl_x_hex_27;
  delete[] ctrl_y_hex_27;
  delete[] ctrl_z_hex_27;
  delete[] ctrl_x_quad_3D;
  delete[] ctrl_y_quad_3D;
  delete[] ctrl_z_quad_3D;
  delete[] sum2; delete[] sum3;
  delete[] sum_quad_3D;
  delete[] R_matlab;
  delete[] dR_dx_matlab; 
  delete[] dR_dy_matlab; 
  delete[] dR_dz_matlab; 
  delete[] d2R_dxx_matlab;
  delete[] d2R_dyy_matlab;
  delete[] d2R_dzz_matlab;
  delete[] d2R_dxy_matlab;
  delete[] d2R_dxz_matlab;
  delete[] d2R_dyz_matlab;
  delete[] J_matlab;
  delete[] Jinv_matlab;
  delete[] basis;
  delete[] dR_dx;
  delete[] dR_dy;
  delete[] dR_dz;
  delete[] d2R_dxx;
  delete[] d2R_dyy;
  delete[] d2R_dzz;
  delete[] d2R_dxy;
  delete[] d2R_dxz;
  delete[] d2R_dyz;
  delete[] J;
  delete[] Jinv;
  delete[] ctrl_x_quad_9;
  delete[] ctrl_y_quad_9;
  delete[] R_quad_matlab;
  delete[] dR_dx_quad_matlab;
  delete[] dR_dy_quad_matlab;
  delete[] d2R_dxx_quad_matlab;
  delete[] d2R_dyy_quad_matlab;
  delete[] d2R_dxy_quad_matlab;
  delete[] J_quad_matlab;
  delete[] Jinv_quad_matlab;
  delete[] ctrl_x_quad9_3D;
  delete[] ctrl_y_quad9_3D;
  delete[] R_quad_3D_matlab;
  delete[] outn_matlab;
  return EXIT_SUCCESS;
}
