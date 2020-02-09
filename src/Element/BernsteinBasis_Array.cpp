#include "BernsteinBasis_Array.hpp"

BernsteinBasis_Array::BernsteinBasis_Array( const int &input_degree,
    const IQuadPts * const &in_quaInfo )
: degree(input_degree), nQuapts(in_quaInfo->get_num_quadPts()),
  degp1(degree+1)
{
  der0 = new double [degp1*nQuapts];
  der1 = new double [degp1*nQuapts];
  der2 = new double [degp1*nQuapts];
  
  for(int ii=0; ii<nQuapts; ++ii)
  {
    const double qp = in_quaInfo->get_qp(ii);
    funcInit_switch(ii, qp);
  }
}


BernsteinBasis_Array::BernsteinBasis_Array( const int &input_degree,
            const double &in_qp )
: degree(input_degree), nQuapts(1), degp1(degree+1)
{
  der0 = new double [degp1];
  der1 = new double [degp1];
  der2 = new double [degp1];
  
  const double qp = in_qp;
  funcInit_switch(0, qp);
}


BernsteinBasis_Array::~BernsteinBasis_Array()
{
  delete [] der0; delete [] der1; delete [] der2;
  der0 = NULL; der1 = NULL; der2 = NULL;
}


void BernsteinBasis_Array::funcInit_switch(const int &quaindex, const double &xi)
{
  const int pos = quaindex * degp1;
  switch(degree)
  {
    case 0:
      der0[pos] = 1.0;
      der1[pos] = 0.0;
      der2[pos] = 0.0;
      break;

    case 1:
      der0[pos+0] = 1.0 - xi;
      der0[pos+1] = xi;
      
      der1[pos+0] = -1.0;
      der1[pos+1] = 1.0;
      
      der2[pos+0] = 0.0;
      der2[pos+1] = 0.0;
      break;

    case 2:
      der0[pos+0] = (1.0-xi) * (1.0-xi);
      der0[pos+1] = 2 * xi * (1.0-xi);
      der0[pos+2] = xi * xi;
      
      der1[pos+0] = 2.0 * (xi - 1);
      der1[pos+1] = 2.0 - 4.0 * xi;
      der1[pos+2] = 2 * xi;
      
      der2[pos+0] = 2.0;
      der2[pos+1] = -4.0;
      der2[pos+2] = 2.0;
      break;
    
    case 3:
      der0[pos+0] = (1.0 - xi) * (1.0 - xi) * (1.0 - xi);
      der0[pos+1] = 3.0 * xi * (1.0 - xi) * (1.0 - xi);
      der0[pos+2] = 3.0 * xi * xi * (1.0 - xi);
      der0[pos+3] = xi * xi * xi;

      der1[pos+0] = -3.0 * (1.0 - xi) * (1.0 - xi);
      der1[pos+1] = 9.0 * xi * xi - 12.0 * xi + 3.0;
      der1[pos+2] = 6.0 * xi - 9.0 * xi * xi;
      der1[pos+3] = 3.0 * xi * xi;

      der2[pos+0] = 6.0 - 6.0 * xi;
      der2[pos+1] = 18.0 * xi - 12.0;
      der2[pos+2] = 6.0 - 18.0 * xi;
      der2[pos+3] = 6.0 * xi;
      break;

    default:
      // ii = 0
      der0[pos] = pow(1.0-xi, degree);
      der1[pos] = -1.0 * degree * pow(1.0-xi, degree-1);
      der2[pos] = degree * (degree-1) * pow(1.0-xi, degree-2);
      
      // ii = 1
      der0[pos+1] = degree * xi * pow(1.0-xi, degree-1);
      der1[pos+1] = degree * pow(1.0-xi, degree-2) * (1.0 - degree * xi);
      der2[pos+1] = degree * pow(1.0-xi, degree-3) * (degree-1) * (degree*xi-2.0);

      for(int ii=2; ii<degree-1; ++ii)
      {
        const double bicoeff = MATH_T::binomialCoefficient(degree, ii);
        der0[pos+ii] = bicoeff * pow(xi, ii) * pow(1.0-xi, degree - ii);
        der1[pos+ii] = bicoeff * pow(xi, ii-1) * pow(1.0 - xi, degree-ii-1) * (ii - degree*xi);
        der2[pos+ii] = bicoeff * pow(xi, ii-2) * pow(1.0-xi, degree-ii-2) * ((degree*degree-degree)*xi*xi + 2*ii*(1.0-degree)*xi + ii*(ii-1.0));
      }
      
      // ii = degree -1
      der0[pos+degree-1] = degree * pow(xi, degree-1) * (1.0-xi);
      der1[pos+degree-1] = degree * pow(xi, degree-2) * (degree-1.0 - degree * xi);
      der2[pos+degree-1] = degree * (degree-1) * pow(xi, degree-3) * (degree-2.0 - degree * xi);

      // ii = degree
      der0[pos+degree] = pow(xi, degree);
      der1[pos+degree] = degree * pow(xi, degree-1);
      der2[pos+degree] = degree * (degree-1) * pow(xi, degree-2);
      break;
  }
}


void BernsteinBasis_Array::print_info() const
{
  std::cout<<"Bernstein Basis: \n";
  std::cout<<"  degree = "<<degree<<'\t';
  std::cout<<"  nQuapts = "<<nQuapts<<'\t';
  std::cout<<"  degp1 = "<<degp1<<'\n';
  for(int ii=0; ii<degp1*nQuapts; ++ii)
    std::cout<<der0[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<degp1*nQuapts; ++ii)
    std::cout<<der1[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<degp1*nQuapts; ++ii)
    std::cout<<der2[ii]<<'\t';
  std::cout<<std::endl;
}

// EOF
