// This is a simple driver for testing the output of the membrane element
// routines.
#include "QuadPts_debug.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"

void print_2Darray(const double * const arr, const int nrow,
  const int ncol);

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const int dim   = 3;
  const int numpt = 4;

  std::vector<double> in_qp{ 1.0/3.0, 1.0/3.0, 1.0/3.0,
                                 0.6,     0.2,     0.2,
                                 0.2,     0.6,     0.2,
                                 0.2,     0.2,     0.6 };
  std::vector<double> in_qw{ -0.5625,
                             0.520833333333333,
                             0.520833333333333,
                             0.520833333333333 };
  
  IQuadPts * quad = new QuadPts_debug(dim, numpt, in_qp, in_qw);

  quad -> print_info();

  FEAElement * elem = new FEAElement_Triangle3_membrane( numpt );
  const int nLocBas = elem -> get_nLocBas();

  double ctrl_x[3] = {  0.3971,  0.4969, 0.4516 };
  double ctrl_y[3] = { -1.4233, -1.2942, 1.3001 };
  double ctrl_z[3] = {  9.7337,  9.6558, 9.8612 };

  elem -> buildBasis( quad, ctrl_x, ctrl_y, ctrl_z );

  // Basis function gradients with respect to lamina coords
  double dR_dxl [nLocBas] = {0.0};
  double dR_dyl [nLocBas] = {0.0};
  
  elem -> get_gradR(0, dR_dxl, dR_dyl);

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
  for(int ii = 0; ii < nLocBas*dim; ++ii)
  {
    for(int jj = 0; jj < nLocBas*dim; ++jj)
    {
      for(int kk = 0; kk < 5; ++kk)
      {
        for(int ll = 0; ll < 5; ++ll)
        {
          Kl[ii*(nLocBas*dim) + jj] +=
            Bl[kk*(nLocBas*dim)+ii] * D[5*kk+ll] * Bl[ll*(nLocBas*dim)+jj];
        }
      }
    }
  }
  std::cout << "\n===== K in lamina coords =====" << std::endl;
  print_2Darray(Kl, nLocBas*dim, nLocBas*dim);

  // Global-to-local rotation matrix Q
  Matrix_3x3 Q = Matrix_3x3();
  elem -> get_rotationMatrix(0, Q);
  std::cout << "\n===== Q =====" << std::endl;
  Q.print();
 
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
          for(int kk = 0; kk < dim; ++kk)
          {
            for(int ll = 0; ll < dim; ++ll)
            {
              Kg[(A*dim+ii)*(nLocBas*dim) + (B*dim+jj)] +=
                Q(kk,ii) * Kl[(A*dim+kk)*(nLocBas*dim) + (B*dim+ll)] * Q(ll, jj);
            }
          }
        }
      }
    }
  }
  std::cout << "\n===== K in global coords =====" << std::endl;
  print_2Darray(Kg, nLocBas*dim, nLocBas*dim);

  delete quad;
  delete elem;

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
      std::cout << std::setw(12) << arr[ii * ncol + jj] << "\t";
    }
    std::cout << std::endl;
  }
}
// EOF
