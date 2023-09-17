#include "Math_Tools.hpp"

double MATH_T::L2Proj_DGP0( const double * const &f,
    const double * const &gwts, const int &nqp )
{
  double sum_top = 0.0, sum_bot = 0.0;
  for(int ii=0; ii<nqp; ++ii)
  {
    sum_top += f[ii] * gwts[ii];
    sum_bot += gwts[ii];
  }
  return sum_top / sum_bot;
}

void MATH_T::L2Proj_DGP1_2D( const double * const &f,
    const double * const &gwts,
    const double * const &qp_x,
    const double * const &qp_y,
    const int &nqp,
    double &coeff_0, double &coeff_x, double &coeff_y )
{
  double a [9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double b [3] = {0.0, 0.0, 0.0};

  for(int ii=0; ii<nqp; ++ii)
  {
    a[0] += gwts[ii];
    a[1] += gwts[ii] * qp_x[ii];
    a[2] += gwts[ii] * qp_y[ii];
    a[4] += gwts[ii] * qp_x[ii] * qp_x[ii];
    a[5] += gwts[ii] * qp_x[ii] * qp_y[ii];
    a[8] += gwts[ii] * qp_y[ii] * qp_y[ii];

    b[0] += gwts[ii] * f[ii];
    b[1] += gwts[ii] * f[ii] * qp_x[ii];
    b[2] += gwts[ii] * f[ii] * qp_y[ii];
  }
  a[3] = a[1]; a[6] = a[2]; a[7] = a[5];

  Matrix_double_3by3_Array AA(a[0], a[1], a[2], a[3],
      a[4], a[5], a[6], a[7], a[8]);

  AA.LU_fac();

  AA.LU_solve(b[0], b[1], b[2], coeff_0, coeff_x, coeff_y);
}

void MATH_T::L2Proj_DGP1_3D( const double * const &f,
    const double * const &gwts,
    const double * const &qp_x,
    const double * const &qp_y,
    const double * const &qp_z,
    const int &nqp,
    double &coeff_0, double &coeff_x, double &coeff_y, double &coeff_z )
{
  double a [16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double b [4]  = {0.0, 0.0, 0.0, 0.0};

  for(int ii=0; ii<nqp; ++ii)
  {
    a[0] += gwts[ii];
    a[1] += gwts[ii] * qp_x[ii];
    a[2] += gwts[ii] * qp_y[ii];
    a[3] += gwts[ii] * qp_z[ii];
    a[5] += gwts[ii] * qp_x[ii] * qp_x[ii];
    a[6] += gwts[ii] * qp_x[ii] * qp_y[ii];
    a[7] += gwts[ii] * qp_x[ii] * qp_z[ii];
    a[10] += gwts[ii] * qp_y[ii] * qp_y[ii];
    a[11] += gwts[ii] * qp_y[ii] * qp_z[ii];
    a[15] += gwts[ii] * qp_z[ii] * qp_z[ii];

    b[0] += gwts[ii] * f[ii];
    b[1] += gwts[ii] * f[ii] * qp_x[ii];
    b[2] += gwts[ii] * f[ii] * qp_y[ii];
    b[3] += gwts[ii] * f[ii] * qp_z[ii];
  }
  a[4] = a[1]; a[8] = a[2]; a[12] = a[3];
  a[9] = a[6]; a[13] = a[7]; a[14] = a[11];

  MATH_T::Matrix_dense AA(4);

  AA.set_values( a );

  AA.LU_fac();

  double out[4] = {0.0, 0.0, 0.0, 0.0};

  AA.LU_solve(b, out);

  coeff_0 = out[0]; coeff_x = out[1]; coeff_y = out[2]; coeff_z = out[3];
}

void MATH_T::get_n_from_t( 
    const double &tx, const double &ty, const double &tz,
    const double &p0_x, const double &p0_y, const double &p0_z,
    const double &p1_x, const double &p1_y, const double &p1_z,
    double &nx, double &ny, double &nz )
{
  const double mx = p0_x - p1_x;
  const double my = p0_y - p1_y;
  const double mz = p0_z - p1_z;

  const double mdt = mx * tx + my * ty + mz * tz;
  const double tdt = tx * tx + ty * ty + tz * tz;
  const double fac = mdt / tdt;

  nx = mx - fac * tx;
  ny = my - fac * ty;
  nz = mz - fac * tz;

  const double len = std::sqrt(nx*nx + ny*ny + nz*nz);
  nx = nx / len;
  ny = ny / len;
  nz = nz / len;
}

void MATH_T::get_tet_sphere_info( const double &x0, const double &x1,
    const double &x2, const double &x3, const double &y0,
    const double &y1, const double &y2, const double &y3,
    const double &z0, const double &z1, const double &z2,
    const double &z3, double &x, double &y, double &z, double &r )
{
  Matrix_double_3by3_Array AA(
      2.0 * (x1-x0), 2.0 * (y1-y0), 2.0 * (z1-z0),
      2.0 * (x2-x0), 2.0 * (y2-y0), 2.0 * (z2-z0),
      2.0 * (x3-x0), 2.0 * (y3-y0), 2.0 * (z3-z0) );

  AA.LU_fac();

  const double xyz2 = x0*x0 + y0*y0 + z0*z0;

  AA.LU_solve( x1*x1 + y1*y1 + z1*z1 - xyz2,
      x2*x2 + y2*y2 + z2*z2 - xyz2,
      x3*x3 + y3*y3 + z3*z3 - xyz2,
      x, y, z );

  r = std::sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) );
}

Vector_3 MATH_T::get_tet_sphere_info( const Vector_3 &pt0,
    const Vector_3 &pt1, const Vector_3 &pt2, const Vector_3 &pt3, 
    double &radius ) 
{
  Matrix_double_3by3_Array AA(
      2.0 * (pt1.x()-pt0.x()), 2.0 * (pt1.y()-pt0.y()), 2.0 * (pt1.z()-pt0.z()),
      2.0 * (pt2.x()-pt0.x()), 2.0 * (pt2.y()-pt0.y()), 2.0 * (pt2.z()-pt0.z()),
      2.0 * (pt3.x()-pt0.x()), 2.0 * (pt3.y()-pt0.y()), 2.0 * (pt3.z()-pt0.z()) );

  AA.LU_fac();

  const double xyz2 = pt0.dot_product( pt0 );

  const Vector_3 centre = AA.LU_solve( Vector_3( pt1.dot_product(pt1) - xyz2,
      pt2.dot_product(pt2) - xyz2, pt3.dot_product(pt3) - xyz2 ) );

  radius = ( centre - pt0 ).norm2();

  return centre;
}

double MATH_T::get_std_dev( const std::vector<double> &vec )
{
  double mean_val = MATH_T::get_mean(vec);
  const unsigned int len = vec.size();

  double sum = 0.0; double nn = 0.0;
  for(unsigned int ii=0; ii<len; ++ii)
  {
    sum += (vec[ii] - mean_val) * (vec[ii] - mean_val);
    nn  += 1.0;
  }

  return std::sqrt( sum / nn );
}

void MATH_T::gen_Gaussian( const int &n, const double &mean, 
    const double &std, std::vector<double> &val )
{
  const int m = n + n % 2;
  val.resize(m);
  for ( int ii = 0; ii < m; ii += 2 )
  {
    double x,y,rsq,f;
    do {
      x = 2.0 * rand() / (double)RAND_MAX - 1.0;
      y = 2.0 * rand() / (double)RAND_MAX - 1.0;
      rsq = x * x + y * y;
    }while( rsq >= 1. || rsq == 0. );
    f = std::sqrt( -2.0 * log(rsq) / rsq );
    val[ii]   = mean + std * x * f;
    val[ii+1] = mean + std * y * f;
  }
}

void MATH_T::print_Histogram( const std::vector<double> &val )
{
  const int width = 50;
  int max = 0;

  const double mean = MATH_T::get_mean(val);
  const double std  = MATH_T::get_std_dev(val);

  const double low   = mean - 3.05 * std;
  const double high  = mean + 3.05 * std;
  const double delta =  0.1 * std;

  const int n = (int)val.size();

  const int nbins = (int)((high - low) / delta);
  int* bins = (int*)calloc(nbins,sizeof(int));
  if ( bins != NULL )
  {
    for ( int i = 0; i < n; i++ )
    {
      int j = (int)( (val[i] - low) / delta );
      if ( 0 <= j  &&  j < nbins ) bins[j]++;
    }

    for ( int j = 0; j < nbins; j++ )
      if ( max < bins[j] ) max = bins[j];

    for ( int j = 0; j < nbins; j++ )
    {
      printf("(%5.2f, %5.2f) |", low + j * delta, low + (j + 1) * delta );
      int k = (int)( (double)width * (double)bins[j] / (double)max );
      while(k-- > 0) putchar('*');
      printf("  %-.1f%%", bins[j] * 100.0 / (double)n);
      putchar('\n');
    }
    free(bins);
  }
}

namespace MATH_T
{
  Matrix_dense::Matrix_dense()
    : N(0), NN(N*N)
  {
    mat = NULL;
    p = NULL;
    invm = NULL;
    is_fac = false;
  }

  Matrix_dense::Matrix_dense(const int &in_msize)
    : N(in_msize), NN(N*N)
  {
    mat = new double [NN];
    p = new int [N];
    invm = new double [N];

    for(int ii=0; ii<NN; ++ii)
      mat[ii] = 0.0;

    for(int ii=0; ii<N; ++ii)
    {
      mat[ii*N+ii] = 1.0;
      p[ii] = ii;
      invm[ii] = 1.0; 
    }

    is_fac = false;
  }

  Matrix_dense::Matrix_dense( const Matrix_dense * const &inm )
    : N(inm->get_size()), NN(N*N)
  {
    mat = new double [NN];
    p = new int [N];
    invm = new double [N];

    for(int ii=0; ii<N; ++ii)
    {
      for(int jj=0; jj<N; ++jj)
        mat[ii*N+jj] = inm->get_mat(ii, jj);

      p[ii] = inm->get_p(ii);
      invm[ii] = inm->get_invm(ii);
    }

    is_fac = false;
  }

  Matrix_dense::~Matrix_dense()
  {
    if(mat != NULL)
    {
      delete [] mat;
      mat = NULL;
    }
    if(p != NULL)
    {
      delete [] p;
      p = NULL;
    }
    if(invm != NULL)
    {
      delete [] invm;
      invm = NULL;
    }
  }

  void Matrix_dense::print_mat(const int &pre) const
  {
    std::cout<<std::setprecision(pre);
    int counter = -1; 
    for(int ii=0; ii<N; ++ii)
    {
      for(int jj=0; jj<N; ++jj)
      {
        counter += 1;
        std::cout<<mat[counter]<<'\t';
      }
      std::cout<<'\n';
    }
  }

  void Matrix_dense::print_info() const
  {
    std::cout<<"N = "<<N<<std::endl;
    std::cout<<"Matrix : "<<std::endl;
    print_mat();
    std::cout<<"p : \n";
    for(int ii=0; ii<N; ++ii)
      std::cout<<p[ii]<<'\t';
    std::cout<<std::endl;

    std::cout<<"invm : \n";
    for(int ii=0; ii<N; ++ii)
      std::cout<<invm[ii]<<'\t';
    std::cout<<std::endl;
    if(is_fac)
      std::cout<<"Matrix is factorized.\n";
    else
      std::cout<<"Matrix is NOT factorized.\n";
  }

  int Matrix_dense::get_size() const
  {
    return N;
  }

  double Matrix_dense::get_mat(const int &ii, const int &jj) const
  {
    return mat[ii*N+jj];
  }

  int Matrix_dense::get_p( const int &ii ) const
  {
    return p[ii];
  }

  double Matrix_dense::get_invm(const int &ii) const
  {
    return invm[ii];
  }

  void Matrix_dense::zero_mat()
  {
    for(int ii=0; ii<NN; ++ii)
      mat[ii] = 0.0;
  }

  void Matrix_dense::set_value(const int &ii, const int &jj, const double &val)
  {
    mat[ii*N+jj] = val;
  }

  void Matrix_dense::set_values( int const * const &index_i, int const * const &index_j,
      double const * const &vals, const int &num )
  {
    for(int ii=0; ii<num; ++ii)
    {
      mat[index_i[ii] * N + index_j[ii]] = vals[ii];
    }
  }

  void Matrix_dense::set_values( double const * const &vals )
  {
    for(int ii=0; ii<NN; ++ii)
      mat[ii] = vals[ii];
  }

  void Matrix_dense::gen_rand()
  {
    srand(time(NULL));

    for(int ii=0; ii<NN; ++ii)
    {
      double value = rand() % 1000; 

      mat[ii] = value * 1.0e-2 - 5;

      p[ii] = ii;
    }
  }

  void Matrix_dense::gen_hilb()
  {
    for(int ii=0; ii<N; ++ii)
    {
      for(int jj=0; jj<N; ++jj)
      {
        mat[ii*N+jj] = 1.0 / (ii+jj+1.0);
      }
      p[ii] = ii;
    }
  }

  void Matrix_dense::gen_id()
  {
    for(int ii=0; ii<NN; ++ii) mat[ii] = 0.0;

    for(int ii=0; ii<N; ++ii){
      mat[ii*N+ii] = 1.0;
      p[ii] = ii;
    }
  }

  void Matrix_dense::LU_fac()
  {
    double max_value, temp, invAkk;
    int max_index, int_temp;
    bool pivot_flag;

    for(int kk=0; kk<N-1; ++kk)
    {
      max_value = std::abs(mat[kk*N+kk]);
      max_index = kk;
      pivot_flag = false;
      for(int ii=kk+1; ii<N; ++ii)
      {
        if( max_value < std::abs(mat[ii*N+kk]) )
        {
          max_value = std::abs(mat[ii*N+kk]);
          max_index = ii;
          pivot_flag = true;
        }
      }

      if(pivot_flag)
      {
        int_temp = p[kk];
        p[kk] = p[max_index];
        p[max_index] = int_temp;

        for(int ii=0; ii<N; ++ii)
        {
          temp = mat[kk*N+ii];
          mat[kk*N+ii] = mat[max_index*N+ii];
          mat[max_index*N+ii] = temp;
        }
      }

      invAkk = 1.0 / mat[kk*N+kk];

      for(int ii=kk+1; ii<N; ++ii)
      {
        mat[ii*N+kk] = mat[ii*N+kk] * invAkk;
        for(int jj=kk+1; jj<N; ++jj)
          mat[ii*N+jj] -= mat[ii*N+kk] * mat[kk*N+jj];
      }
    }

    is_fac = true;

    for(int kk=0; kk<N; ++kk)
      invm[kk] = 1.0 / mat[kk*N+kk];
  }

  void Matrix_dense::LU_solve( double const * const &b, double * const &x ) const
  {
    for(int ii=0; ii<N; ++ii)
      x[ii] = b[p[ii]];

    for(int ii=1; ii<N; ++ii)
      for(int jj=0; jj<ii; ++jj)
        x[ii] -= mat[ii*N+jj] * x[jj];

    for(int ii=N-1; ii>=0; --ii)
    {
      for(int jj=N-1; jj>ii; --jj)
        x[ii] -= mat[ii*N+jj] * x[jj];
      x[ii] = x[ii] * invm[ii];
    }
  }

  void Matrix_dense::Axb( double const * const &b, double * const &x ) const
  {
    if( is_fac ) std::cout<<"Warning: this matrix has been factorized.\n";

    for(int ii=0; ii<N; ++ii)
    {
      x[ii] = 0.0;
      for(int jj=0; jj<N; ++jj)
        x[ii] += mat[ii*N+jj] * b[jj];
    }
  }


  Matrix_SymPos_dense::Matrix_SymPos_dense() : Matrix_dense()
  {}


  Matrix_SymPos_dense::Matrix_SymPos_dense( const int &msize )
    : Matrix_dense(msize)
  {}


  Matrix_SymPos_dense::Matrix_SymPos_dense( 
      const Matrix_SymPos_dense * const &in_mat )
    : Matrix_dense(in_mat)
  {}


  Matrix_SymPos_dense::~Matrix_SymPos_dense()
  {}

  void Matrix_SymPos_dense::check_symm() const
  {
    for(int ii = 0; ii<N; ++ii)
    {
      for(int jj = 0; jj<ii; ++jj)
      {
        if( !MATH_T::equals( mat[ii*N+jj], mat[jj*N+ii], 1.0e-15) ) 
          std::cout<<"error: Matrix_SymPos entry ("<<ii<<","<<jj<<") does not match entry (jj"<<","<<ii<<"). \n";
      }
    }
  }

  void Matrix_SymPos_dense::copy( const Matrix_SymPos_dense * const &in_mat )
  {
    for(int ii=0; ii<NN; ++ii) mat[ii] = in_mat->get_mat(ii);

    for(int ii=0; ii<N; ++ii)
    {
      p[ii] = in_mat->get_p(ii);
      invm[ii] = in_mat->get_invm(ii);
    }
    is_fac = false; 
  }

  void Matrix_SymPos_dense::LDLt_fac()
  {
    // This algorithm is given in Shufang XU's book, pp 31. 
    double v[N]; 

    for(int jj=0; jj<N; ++jj)
    {
      const int Njj = jj * N;
      for(int kk=0; kk<jj; ++kk) v[kk] = mat[Njj+kk] * mat[kk*N+kk];

      for(int kk=0; kk<jj; ++kk) mat[Njj+jj] -= v[kk] * mat[Njj+kk];

      invm[jj] = 1.0 / mat[Njj+jj];

      for(int ii=jj+1; ii<N; ++ii)
      {
        for(int kk=0; kk<jj; ++kk) mat[N*ii+jj] -= mat[ii*N+kk] * v[kk];

        mat[N*ii+jj] *= invm[jj];
      }
    }

    is_fac = true;
  }

  void Matrix_SymPos_dense::LDLt_solve( double const * const &b,
      double * const &x ) const
  {
    // Solve for Ly = b
    for(int ii=0; ii<N; ++ii)
    {
      x[ii] = b[ii];
      for(int jj=0; jj<ii; ++jj) x[ii] -= mat[ii*N+jj] * x[jj];
    }

    // Solve for D z = y;
    for(int ii=0; ii<N; ++ii) x[ii] *= invm[ii];

    // Solve L^t x = z
    for(int ii=N-2; ii>=0; --ii)
    {
      for(int jj=ii+1; jj<N; ++jj) x[ii] -= mat[jj*N+ii] * x[jj];
    }
  }
}
// End of Namespace Math_T
// EOF
