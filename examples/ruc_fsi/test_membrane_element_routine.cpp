// This is a simple driver for testing the output of the membrane element
// routines.
#include "QuadPts_debug.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "FEAElement_Triangle6.hpp"

void print_2Darray(const double * const arr, const int nrow,
  const int ncol);

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  const int nLocBas = 3;
  const int dim     = 3;
  int numpt;

  std::vector<double> in_qp;
  std::vector<double> in_qw;

  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  double * sol_wall_disp = new double [dim * nLocBas];

  if(nLocBas == 3)
  {
    numpt = 4;

    in_qp = { 1.0/3.0, 1.0/3.0, 1.0/3.0,
                  0.6,     0.2,     0.2,
                  0.2,     0.6,     0.2,
                  0.2,     0.2,     0.6 };
    in_qw = { -0.5625,
              0.520833333333333,
              0.520833333333333,
              0.520833333333333 };

    ctrl_x[0] =  0.3971; ctrl_x[1] =  0.4969; ctrl_x[2] = 0.4516;
    ctrl_y[0] = -1.4233; ctrl_y[1] = -1.2942; ctrl_y[2] = 1.3001;
    ctrl_z[0] =  9.7337; ctrl_z[1] =  9.6558; ctrl_z[2] = 9.8612;

    sol_wall_disp[0] =  1.0; sol_wall_disp[1] = 2.134; sol_wall_disp[2] = -3.14;
    sol_wall_disp[3] = -1.0; sol_wall_disp[4] = 1.156; sol_wall_disp[5] = -1.56;
    sol_wall_disp[6] =  2.0; sol_wall_disp[7] =  -3.9; sol_wall_disp[8] = 0.271; 
  }
  else if(nLocBas == 6)
  {
    numpt = 13;

    const double a = 0.638444188569809;
    const double b = 0.312865496004875;
    const double c = 1.0 - a - b;
    const double w = 0.077113760890257;
    in_qp = {           1.0/3.0,           1.0/3.0,           1.0/3.0,
              0.479308067841923, 0.260345966079038, 0.260345966079038,
              0.260345966079038, 0.479308067841923, 0.260345966079038,
              0.260345966079038, 0.260345966079038, 0.479308067841923,
              0.869739794195568, 0.065130102902216, 0.065130102902216,
              0.065130102902216, 0.869739794195568, 0.065130102902216,
              0.065130102902216, 0.065130102902216, 0.869739794195568,
                              a,                 b,                 c,
                              a,                 c,                 b,
                              b,                 a,                 c,
                              c,                 a,                 b,
                              b,                 c,                 a,
                              c,                 b,                 a };
    in_qw = {-0.149570044467670,
              0.175615257433204,
              0.175615257433204,
              0.175615257433204,
              0.053347235608839,
              0.053347235608839,
              0.053347235608839,
                              w,
                              w,
                              w,
                              w,
                              w,
                              w };

    // ====== Identical to Triangle3 test geometry with mid-side nodes ======
    // ctrl_x[0] =   0.3971; ctrl_x[1] =  0.4969; ctrl_x[2] =  0.4516;
    // ctrl_x[3] =   0.4470; ctrl_x[4] = 0.47425; ctrl_x[5] = 0.42435;
    // ctrl_y[0] =  -1.4233; ctrl_y[1] = -1.2942; ctrl_y[2] =  1.3001;
    // ctrl_y[3] = -1.35875; ctrl_y[4] = 0.00295; ctrl_y[5] = -0.0616;
    // ctrl_z[0] =   9.7337; ctrl_z[1] =  9.6558; ctrl_z[2] =  9.8612;
    // ctrl_z[3] =  9.69475; ctrl_z[4] =  9.7585; ctrl_z[5] = 9.79745;

    // ====== Curved triangle in xy-plane for comparison with Triangle6 ====== 
    ctrl_x[0] = 8.1622; ctrl_x[1] = 7.8695; ctrl_x[2] = 7.9049;
    ctrl_x[3] = 8.0127; ctrl_x[4] = 7.8685; ctrl_x[5] = 8.0386;
    ctrl_y[0] = 5.0793; ctrl_y[1] = 5.0431; ctrl_y[2] = 5.2765;
    ctrl_y[3] = 5.0724; ctrl_y[4] = 5.1655; ctrl_y[5] = 5.1846;
    ctrl_z[0] =    0.0; ctrl_z[1] =    0.0; ctrl_z[2] =    0.0;
    ctrl_z[3] =    0.0; ctrl_z[4] =    0.0; ctrl_z[5] =    0.0;

    sol_wall_disp[0]  =  1.0; sol_wall_disp[1]  = 2.134; sol_wall_disp[2]  = -3.14;
    sol_wall_disp[3]  = -1.0; sol_wall_disp[4]  = 1.156; sol_wall_disp[5]  = -1.56;
    sol_wall_disp[6]  =  2.0; sol_wall_disp[7]  =  -3.9; sol_wall_disp[8]  = 0.271; 
    sol_wall_disp[9]  =  1.0; sol_wall_disp[10] = 2.134; sol_wall_disp[11] = -3.14;
    sol_wall_disp[12] = -1.0; sol_wall_disp[13] = 1.156; sol_wall_disp[14] = -1.56;
    sol_wall_disp[15] =  2.0; sol_wall_disp[16] =  -3.9; sol_wall_disp[17] = 0.271; 
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");
 
  IQuadPts * quad = new QuadPts_debug(in_qp, in_qw, dim);

  quad -> print_info();

  std::cout << "\n===== sol_wall_disp =====" << std::endl;
  print_2Darray(sol_wall_disp, nLocBas*dim, 1);

  FEAElement * elem = nullptr;
  if(nLocBas == 3) elem = new FEAElement_Triangle3_membrane( numpt );
  else             elem = new FEAElement_Triangle6_membrane( numpt );

  elem -> buildBasis( quad, ctrl_x, ctrl_y, ctrl_z );

  // Basis function gradients with respect to lamina coords
  double dR_dxl [nLocBas] = {0.0};
  double dR_dyl [nLocBas] = {0.0};
  
  // Test on a single quadrature point
  const int qua = 0;

  elem -> get_gradR(qua, dR_dxl, dR_dyl);

  // Global-to-local rotation matrix Q
  const Tensor2_3D Q = elem -> get_rotationMatrix(0);
  std::cout << "\n===== Q =====" << std::endl;
  Q.print();

  double sol_wall_disp_l [ nLocBas * dim ] = {0.0};
  for(int ii=0; ii<nLocBas; ++ii)
  {
    sol_wall_disp_l[dim*ii]   = sol_wall_disp[dim*ii] * Q(0, 0) 
      + sol_wall_disp[dim*ii+1] * Q(0, 1) + sol_wall_disp[dim*ii+2] * Q(0, 2); 

    sol_wall_disp_l[dim*ii+1] = sol_wall_disp[dim*ii] * Q(1, 0) 
      + sol_wall_disp[dim*ii+1] * Q(1, 1) + sol_wall_disp[dim*ii+2] * Q(1, 2); 

    sol_wall_disp_l[dim*ii+2] = sol_wall_disp[dim*ii] * Q(2, 0) 
      + sol_wall_disp[dim*ii+1] * Q(2, 1) + sol_wall_disp[dim*ii+2] * Q(2, 2); 
  }

  std::cout << "\n===== sol_wall_disp in lamina coords =====" << std::endl;
  print_2Darray(sol_wall_disp_l, nLocBas*dim, 1);

  // Strain displacement matrix B in lamina coords
  // 5 x (nLocBas * dim)
  double Bl [5 * nLocBas * dim] = {0.0};
  for(int ii = 0; ii < nLocBas; ++ii)
  {
    Bl[0*nLocBas*dim + ii*dim]     = dR_dxl[ii]; // u1,1
    Bl[1*nLocBas*dim + ii*dim + 1] = dR_dyl[ii]; // u2,2
    Bl[2*nLocBas*dim + ii*dim]     = dR_dyl[ii]; // u1,2
    Bl[2*nLocBas*dim + ii*dim + 1] = dR_dxl[ii]; // u2,1
    Bl[3*nLocBas*dim + ii*dim + 2] = dR_dxl[ii]; // u3,1
    Bl[4*nLocBas*dim + ii*dim + 2] = dR_dyl[ii]; // u3,2
  }
  std::cout << "\n====== B in lamina coords ======" << std::endl;
  print_2Darray(Bl, 5, nLocBas * dim);

  // Elasticity tensor D
  const double nu    = 0.5;
  const double kappa = 0.833333;
  const double E     = 2500000;
  const double coef  = E / (1.0 - nu*nu);

  double D[5 * 5] = {0.0};
  D[0*5 + 0] = coef * 1.0;
  D[0*5 + 1] = coef * nu;
  D[1*5 + 0] = coef * nu;
  D[1*5 + 1] = coef * 1.0;
  D[2*5 + 2] = coef * (1.0 - nu) / 2.0;
  D[3*5 + 3] = coef * kappa * (1.0 - nu) / 2.0;
  D[4*5 + 4] = coef * kappa * (1.0 - nu) / 2.0;
  std::cout << "\n===== D =====" << std::endl;
  print_2Darray(D, 5, 5);

  // Stiffness tensor in lamina coords
  // Bl^T * D * Bl = Bl_{ki} * D_{kl} * Bl_{lj}
  double Kl [(nLocBas*dim) * (nLocBas*dim)] = {0.0};

  // Bl * ul = Bl_{ij} * ul_{j}
  double Blul [5] = {0.0};
  for(int ii=0; ii<5; ++ii)
  {
    for(int jj=0; jj<nLocBas*dim; ++jj)
      Blul[ii] += Bl[ii*(nLocBas*dim)+jj] * sol_wall_disp_l[jj];
  }
  std::cout << "\n====== Bl * ul in lamina coords ======" << std::endl;
  print_2Darray(Blul, 5, 1);


  // D * Bl * ul = D_{ik} * Bl_{kj} * ul_{j}
  double DBlul [5] = {0.0};
  for(int ii=0; ii<5; ++ii)
  {
    for(int jj=0; jj<nLocBas*dim; ++jj)
    {
      for(int kk=0; kk<5; ++kk)
        DBlul[ii] += D[5*ii+kk] * Bl[kk*(nLocBas*dim)+jj] * sol_wall_disp_l[jj];
    }
  }
  std::cout << "\n====== D * Bl * ul in lamina coords ======" << std::endl;
  print_2Darray(DBlul, 5, 1);

  for(int A = 0; A < nLocBas; ++A)
  {
    const double NA_xl = dR_dxl[A], NA_yl = dR_dyl[A];

    for(int B = 0; B < nLocBas; ++B)
    {
      const double NB_xl = dR_dxl[B], NB_yl = dR_dyl[B];

      // Momentum-x with respect to u1, u2 
      Kl[(nLocBas*dim)*(A*dim) + (B*dim)]     += coef * ( NA_xl * NB_xl
          + 0.5*(1.0-nu) * NA_yl * NB_yl );
      Kl[(nLocBas*dim)*(A*dim) + (B*dim+1)]   += coef * ( nu * NA_xl * NB_yl
          + 0.5*(1.0-nu) * NA_yl * NB_xl );

      // Momentum-y with respect to u1, u2 
      Kl[(nLocBas*dim)*(A*dim+1) + (B*dim)]   += coef * ( nu * NA_yl * NB_xl
          + 0.5*(1.0-nu) * NA_xl * NB_yl );
      Kl[(nLocBas*dim)*(A*dim+1) + (B*dim+1)] += coef * ( NA_yl * NB_yl
          + 0.5*(1.0-nu) * NA_xl * NB_xl );

      // Momentum-z with respect to u3 
      Kl[(nLocBas*dim)*(A*dim+2) + (B*dim+2)] += coef * 0.5*kappa*(1.0-nu) * (
          NA_xl * NB_xl + NA_yl * NB_yl );
    } 
  }

  std::cout << "\n===== K in lamina coords =====" << std::endl;
  print_2Darray(Kl, nLocBas*dim, nLocBas*dim);

  // Stiffness tensor in global coords
  // theta^T * Kl * theta, where theta = [Q, 0, 0; 0, Q, 0; 0, 0, Q]
  // or Q^T * Kl_[AB] * Q = Q_{ki} * Kl_[AB]{kl} * Q_{lj}
  double Kg[(nLocBas*dim) * (nLocBas*dim)] = {0.0};
  for(int A = 0; A < nLocBas; ++A)
  {
    for(int B = 0; B < nLocBas; ++B)
    {
      for(int ii = 0; ii < dim; ++ii)
      {
        for(int jj = 0; jj < dim; ++jj)
        {
          // Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(kk,ii) * Kl[ (A*dim+kk)*(snLocBas*dim) + (B*dim+ll) ] * Q(ll, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(nLocBas*dim) + (B*dim+0) ] * Q(0, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(nLocBas*dim) + (B*dim+1) ] * Q(1, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(nLocBas*dim) + (B*dim+2) ] * Q(2, jj);

          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(nLocBas*dim) + (B*dim+0) ] * Q(0, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(nLocBas*dim) + (B*dim+1) ] * Q(1, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(nLocBas*dim) + (B*dim+2) ] * Q(2, jj);
          
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(nLocBas*dim) + (B*dim+0) ] * Q(0, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(nLocBas*dim) + (B*dim+1) ] * Q(1, jj);
          Kg[ (nLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(nLocBas*dim) + (B*dim+2) ] * Q(2, jj);
        }
      }
    }
  }
  std::cout << "\n===== K in global coords =====" << std::endl;
  print_2Darray(Kg, nLocBas*dim, nLocBas*dim);

  // Multiply by displacements in global coords
  // Kg_{ij} * u_{j}
  double lin_elasticity [ nLocBas * dim ] = {0.0};
  for(int ii=0; ii<nLocBas*dim; ++ii)
    for(int jj=0; jj<nLocBas*dim; ++jj)
      lin_elasticity[ii] += Kg[ (nLocBas*dim)*ii + jj ] * sol_wall_disp[jj];

  std::cout << "\n===== lin_elasticity v1 =====" << std::endl;
  print_2Darray(lin_elasticity, nLocBas*dim, 1);

  // Generate lin_elasticity without first generating the tangent ===========
  double u1l_xl = 0.0, u1l_yl = 0.0;
  double u2l_xl = 0.0, u2l_yl = 0.0;
  double u3l_xl = 0.0, u3l_yl = 0.0;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    u1l_xl += sol_wall_disp_l[dim*ii]   * dR_dxl[ii];
    u1l_yl += sol_wall_disp_l[dim*ii]   * dR_dyl[ii];
    u2l_xl += sol_wall_disp_l[dim*ii+1] * dR_dxl[ii];
    u2l_yl += sol_wall_disp_l[dim*ii+1] * dR_dyl[ii];
    u3l_xl += sol_wall_disp_l[dim*ii+2] * dR_dxl[ii];
    u3l_yl += sol_wall_disp_l[dim*ii+2] * dR_dyl[ii];
  }

  std::cout << "\n====== basis fcn gradients in lamina coords V2 ======" << std::endl;
  std::cout << u1l_xl << " " << u1l_yl << " " << u2l_xl << " " << u2l_yl << " " << u3l_xl << " " << u3l_yl << std::endl;

  double Blul2 [5] = {u1l_xl, u2l_yl, u1l_yl + u2l_xl, u3l_xl, u3l_yl};
  std::cout << "\n====== Bl * ul in lamina coords V2 ======" << std::endl;
  print_2Darray(Blul2, 5, 1);

  // Cauchy stress in lamina coords 
  // double sigma_l [ dim * dim ] = {0.0};
  // sigma_l[0] = coef * (u1l_xl + nu*u2l_yl);              // sigma_l_11
  // sigma_l[1] = coef * 0.5*(1.0-nu) * (u1l_yl + u2l_xl);  // sigma_l_12
  // sigma_l[3] = coef * 0.5*(1.0-nu) * (u1l_yl + u2l_xl);  // sigma_l_21
  // sigma_l[2] = coef * 0.5*kappa*(1.0-nu) * u3l_xl;       // sigma_l_13
  // sigma_l[6] = coef * 0.5*kappa*(1.0-nu) * u3l_xl;       // sigma_l_31  
  // sigma_l[4] = coef * (nu*u1l_xl + u2l_yl);              // sigma_l_22
  // sigma_l[5] = coef * 0.5*kappa*(1.0-nu) * u3l_yl;       // sigma_l_23
  // sigma_l[7] = coef * 0.5*kappa*(1.0-nu) * u3l_yl;       // sigma_l_32 
  // sigma_l[8] = 0.0;                                      // sigma_l_33

  Tensor2_3D sigma = Tensor2_3D(
      u1l_xl + nu*u2l_yl,               0.5*(1.0-nu) * (u1l_yl + u2l_xl), 0.5*kappa*(1.0-nu) * u3l_xl,
      0.5*(1.0-nu) * (u1l_yl + u2l_xl), nu*u1l_xl + u2l_yl,               0.5*kappa*(1.0-nu) * u3l_yl,
      0.5*kappa*(1.0-nu) * u3l_xl,      0.5*kappa*(1.0-nu) * u3l_yl,      0.0 );
  sigma *= coef;

  std::cout << "\n===== sigma_l =====" << std::endl;
  sigma.print();

  // Cauchy stress in global coords
  // Q^T * sigma_l * Q = Q_{ki} * sigma_l_{kl} * Q_{lj}
  sigma.MatRot(Q);

  std::cout << "\n===== sigma_g =====" << std::endl;
  sigma.print();

  // Basis function gradients with respect to global coords
  // dR/dx_{i} = Q_{ji} * dR/dxl_{j}
  // Note that dR/dzl = 0.0
  double dR_dx [nLocBas] = {0.0};
  double dR_dy [nLocBas] = {0.0};
  double dR_dz [nLocBas] = {0.0};

  for(int ii=0; ii<nLocBas; ++ii)
  {
    dR_dx[ii] = Q(0, 0) * dR_dxl[ii] + Q(1, 0) * dR_dyl[ii];
    dR_dy[ii] = Q(0, 1) * dR_dxl[ii] + Q(1, 1) * dR_dyl[ii];
    dR_dz[ii] = Q(0, 2) * dR_dxl[ii] + Q(1, 2) * dR_dyl[ii];
  }
  
  std::cout << "\n===== Triangle6_membrane dR_dxyz =====" << std::endl;
  print_2Darray(dR_dx, nLocBas, 1);
  print_2Darray(dR_dy, nLocBas, 1);
  print_2Darray(dR_dz, nLocBas, 1);

  double lin_elasticity2 [ nLocBas * dim ] = {0.0};
  for(int A=0; A<nLocBas; ++A)
  {
    const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

    lin_elasticity2[dim*A]   = NA_x * sigma(0, 0)
      + NA_y * sigma(0, 1) + NA_z * sigma(0, 2);

    lin_elasticity2[dim*A+1] = NA_x * sigma(1, 0)
      + NA_y * sigma(1, 1) + NA_z * sigma(1, 2);

    lin_elasticity2[dim*A+2] = NA_x * sigma(2, 0)
      + NA_y * sigma(2, 1) + NA_z * sigma(2, 2);
  }

  std::cout << "\n===== lin_elasticity v2 =====" << std::endl;
  print_2Darray(lin_elasticity2, nLocBas*dim, 1);

  // Use Triangle6 basis function gradients to verify Triangle6_membrane
  FEAElement * elem_tri6 = nullptr; 
  if(nLocBas == 6)
  {
    elem_tri6 = new FEAElement_Triangle6( numpt );
    elem_tri6 -> buildBasis( quad, ctrl_x, ctrl_y );

    double dR_dx [nLocBas] = {0.0};
    double dR_dy [nLocBas] = {0.0};

    // Test on a single quadrature point
    const int qua = 0;

    elem_tri6 -> get_gradR(qua, dR_dx, dR_dy);

    std::cout << "\n===== Triangle6 dR_dx =====" << std::endl;
    print_2Darray(dR_dx, nLocBas, 1);

    // Stiffness tensor
    // B^T * D * B = B_{ki} * D_{kl} * B_{lj}
    double K [(nLocBas*dim) * (nLocBas*dim)] = {0.0};
    for(int A = 0; A < nLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A];

      for(int B = 0; B < nLocBas; ++B)
      {
        const double NB_x = dR_dx[B], NB_y = dR_dy[B];

        // Momentum-x with respect to u1, u2 
        K[(nLocBas*dim)*(A*dim) + (B*dim)]     += coef * ( NA_x * NB_x
            + 0.5*(1.0-nu) * NA_y * NB_y );
        K[(nLocBas*dim)*(A*dim) + (B*dim+1)]   += coef * ( nu * NA_x * NB_y
            + 0.5*(1.0-nu) * NA_y * NB_x );

        // Momentum-y with respect to u1, u2 
        K[(nLocBas*dim)*(A*dim+1) + (B*dim)]   += coef * ( nu * NA_y * NB_x
            + 0.5*(1.0-nu) * NA_x * NB_y );
        K[(nLocBas*dim)*(A*dim+1) + (B*dim+1)] += coef * ( NA_y * NB_y
            + 0.5*(1.0-nu) * NA_x * NB_x );

        // Momentum-z with respect to u3 
        K[(nLocBas*dim)*(A*dim+2) + (B*dim+2)] += coef*0.5*kappa*(1.0-nu) * (
            NA_x * NB_x + NA_y * NB_y );
      } 
    }

    std::cout << "\n===== K for Triangle6 =====" << std::endl;
    print_2Darray(K, nLocBas*dim, nLocBas*dim);
  }


  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] sol_wall_disp;
  ctrl_x = nullptr; ctrl_y = nullptr; ctrl_z = nullptr; sol_wall_disp = nullptr;
  delete quad; delete elem;

  if(nLocBas == 6) delete elem_tri6;

  PetscFinalize();
  return EXIT_SUCCESS;
}

void print_2Darray(const double * const arr, const int nrow,
  const int ncol)
{
  for(int ii = 0; ii < nrow; ++ii)
  {
    for(int jj = 0; jj < ncol; ++jj)
    {
      std::cout << std::scientific << std::setprecision(3) << std::setw(10) << arr[ii * ncol + jj] << " ";
      // std::cout << arr[ii * ncol + jj] << " ";
    }
    std::cout << std::endl;
  }
}
// EOF
