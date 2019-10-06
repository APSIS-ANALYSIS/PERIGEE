#include "BernsteinBasis_1D.hpp"

BernsteinBasis_1D::BernsteinBasis_1D( const int in_deg, const IQuadPts * const &in_quaPt )
: deg(in_deg), deg_p1(in_deg + 1), nqua(in_quaPt->get_num_quadPts()),
  vseg(deg_p1 * nqua)
{
  val = new double [3*vseg]; 

  for(int ii=0; ii<nqua; ++ii)
  {
    const double qp = in_quaPt->get_qp(ii);
    eval_val(ii, qp);
  }

}


BernsteinBasis_1D::~BernsteinBasis_1D()
{
  delete [] val; val = NULL;
}


BernsteinBasis_1D::BernsteinBasis_1D()
: deg(0), deg_p1(1), nqua(0), vseg(deg_p1*nqua)
{
  val = NULL;
}


void BernsteinBasis_1D::get_B(const int &qua, std::vector<double> &vec) const
{
  vec.resize(deg_p1);
  const int starter = qua * deg_p1;
  for(int ii=0; ii<deg_p1; ++ii)
    vec[ii] = val[ii + starter];
  VEC_T::shrink2fit(vec);
}


void BernsteinBasis_1D::get_dB_ds(const int &qua, std::vector<double> &vec) const
{
  vec.resize(deg_p1);
  const int starter = qua * deg_p1 + vseg;
  for(int ii=0; ii<deg_p1; ++ii)
    vec[ii] = val[ii + starter];
  VEC_T::shrink2fit(vec);
}


void BernsteinBasis_1D::get_d2B_dss(const int &qua, std::vector<double> &vec) const
{
  vec.resize(deg_p1);
  const int starter = qua * deg_p1 + 2 * vseg;
  for(int ii=0; ii<deg_p1; ++ii)
    vec[ii] = val[ii + starter];
  VEC_T::shrink2fit(vec);
}


void BernsteinBasis_1D::eval_val( const int &quaindex, const double &x )
{
  const int pos = quaindex * deg_p1;
  switch(deg)
  {
    case 1:
      val[pos+0] = 1.0 - x;
      val[pos+1] = x;

      val[vseg + pos + 0] = -1.0;
      val[vseg + pos + 1] = 1.0;

      val[2*vseg + pos + 0 ] = 0.0;
      val[2*vseg + pos + 1 ] = 0.0;
      break;

    case 2:
      val[pos+0] = (1.0-x)*(1.0-x);
      val[pos+1] = 2.0 * x * (1.0 -x);
      val[pos+2] = x*x;

      val[vseg+pos+0] = 2.0 * (x - 1.0);
      val[vseg+pos+1] = 2.0 - 4.0 * x;
      val[vseg+pos+2] = 2.0 * x;

      val[2*vseg+pos+0] = 2.0;
      val[2*vseg+pos+1] = -4.0;
      val[2*vseg+pos+2] = 2.0;
      break;

    case 3:
      val[pos+0] = (1.0-x)*(1.0-x)*(1.0-x);
      val[pos+1] = 3.0 * x * (1.0-x) * (1.0 -x);
      val[pos+2] = 3.0 * x * x * (1.0 -x);
      val[pos+3] = x*x*x;

      val[vseg+pos+0] = -3.0 * (x-1.0) * (x-1.0);
      val[vseg+pos+1] = 9.0 * x * x - 12.0 * x + 3.0;
      val[vseg+pos+2] = -9.0 * x * x + 6.0 * x;
      val[vseg+pos+3] = 3 * x*x;

      val[2*vseg+pos+0] = 6.0 - 6.0 * x;
      val[2*vseg+pos+1] = 18.0*x - 12.0;
      val[2*vseg+pos+2] = -18.0 * x + 6.0;
      val[2*vseg+pos+3] = 6.0 * x;
      break;

    case 4:
      val[pos+0] = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x);
      val[pos+1] = 4.0 * x * (1.0-x) * (1.0-x) * (1.0-x);
      val[pos+2] = 6.0 * x * x * (1.0-x) * (1.0-x);
      val[pos+3] = 4.0 * x * x * x * (1.0 -x);
      val[pos+4] = x * x * x * x; 

      val[vseg+pos+0] = 4.0 * (x-1.0) * (x-1.0) * (x-1.0);
      val[vseg+pos+1] = (1.0-x)*(1.0-x)*(4.0 - 16.0*x);
      val[vseg+pos+2] = 12.0 * x * (x - 1.0) * (2.0 * x - 1.0);
      val[vseg+pos+3] = 12.0 * x * x - 16.0 * x * x * x;
      val[vseg+pos+4] = 4.0 * x * x * x;


      val[2*vseg+pos+0] = 12.0 * (x-1.0) * (x-1.0);
      val[2*vseg+pos+1] = -48.0 * x*x + 72.0 * x - 24.0;
      val[2*vseg+pos+2] = 72.0 * x*x - 72.0 * x + 12.0;
      val[2*vseg+pos+3] = 24.0 * x - 48.0 * x * x;
      val[2*vseg+pos+4] = 12.0 * x * x;
      break;

    case 5:
      val[pos+0] = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(1.0-x);
      val[pos+1] = 5.0 * x * (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x);
      val[pos+2] = 10.0 * x * x * (1.0-x) * (1.0-x) * (1.0-x);
      val[pos+3] = 10.0 * x * x * x * (1.0-x) * (1.0-x);
      val[pos+4] = 5.0 * x * x * x * x * (1.0-x);
      val[pos+5] = x * x * x * x * x;

      val[vseg+pos+0] = -5.0 * (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x);
      val[vseg+pos+1] = 5.0 * (1.0-x)*(1.0-x)*(1.0-x)*(1.0-5.0*x);
      val[vseg+pos+2] = 10.0 * x * (1.0-x)*(1.0-x)*(2.0-5.0*x);
      val[vseg+pos+3] = 10.0 * x * x * (1.0-x)*(3.0-5.0*x);
      val[vseg+pos+4] = 20.0 * x * x * x - 25.0 * x * x * x * x;
      val[vseg+pos+5] = 5.0 * x * x * x * x;

      val[2*vseg+pos+0] = 20.0 * (1.0-x)*(1.0-x)*(1.0-x);
      val[2*vseg+pos+1] = (1.0-x) * (1.0-x) * (100.0*x - 40.0); 
      val[2*vseg+pos+2] = -200.0*x*x*x + 360.0 * x*x - 180.0 * x + 20.0;
      val[2*vseg+pos+3] = 200.0 * x*x*x - 240.0 * x * x + 60.0 * x;
      val[2*vseg+pos+4] = 60.0 * x * x - 100.0 * x * x * x;
      val[2*vseg+pos+5] = 20.0 * x *x * x;
      break;
    
    default:
      SYS_T::commPrint("Error: This Bernstein polynomial is not implemented in BernsteinBasis_1D class. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}

// EOF
