#ifndef MATH_TOOLS_HPP
#define MATH_TOOLS_HPP
// ============================================================================
// Math_Tools.hpp
// This file defines mathematical constants and functions.
//
// Date Created: Feb. 13 2016
// ============================================================================
#include <stdio.h>
#include <vector> 
#include "Matrix_double_3by3_Array.hpp"

namespace MATH_T
{
  // --------------------------------------------------------------------------
  // Useful Constants:
  // --------------------------------------------------------------------------
  // PI = 3.1415926......
  constexpr double PI = 3.14159265358979323846264338327950288419716939937510582;

  // E = 2.71828......
  constexpr double E  = 2.71828182845904523536028747135266249775724709369995957;

  // Compute the binomial coefficient, e.g.
  // n = 0 : 1                  n = 1 : 1 1
  // n = 2 : 1 2 1              n = 3 : 1 3 3 1
  // n = 4 : 1 4 6 4 1          n = 5 : 1 5 10 10 5 1
  inline double binomialCoefficient( const double &n, const double &k )
  {
    if( (k<0) || (k>n) ) return 0.0;

    int m = k;
    if( k > n-k ) m = n-k;

    double c = 1.0;
    for(int ii=1; ii<=m; ++ii) c = c * ( n - m + ii ) / ii;
    return c;
  }

  // ----------------------------------------------------------------
  // Assume ii = iz * dim_x * dim_y + iy * dim_x + ix
  // this function will return ix iy and iz based on the input ii,
  // dim_x, dim_y.
  // ----------------------------------------------------------------
  inline void get_xyz_index( const int &ii, const int &dim_x, const int &dim_y,
      int &ix, int &iy, int &iz)
  {
    const int ixy = ii % (dim_x * dim_y);
    iz = (ii - ixy) / (dim_x * dim_y);
    ix = ixy % dim_x; iy = (ixy - ix) / dim_x;
  }

  // ----------------------------------------------------------------
  // Assume ii = iy * dim_x + ix;
  // this function will return ix and iy based on the input ii and dim_x.
  // ----------------------------------------------------------------
  inline void get_xy_index( const int &ii, const int &dim_x, int &ix, int &iy)
  {
    ix = ii % dim_x;
    iy = (ii-ix)/dim_x;
  }

  // --------------------------------------------------------------------------
  // Useful Functions:
  // --------------------------------------------------------------------------
  // This is the tool to calculate if two double data are within tol
  // default tolerance is 1.0x10^-12. 
  // --------------------------------------------------------------------------
  inline bool equals( const double &a, const double &b, const double &tol=1.0e-12 )
  {
    return ( std::abs(a-b)<tol );
  }
  
  // ----------------------------------------------------------------
  // Generate outward normal vector from a tangential vector.
  // t : the tangential vector
  // p0 : the starting point of the tangential vector
  // p1 : the interior point 
  // n : the normal vector
  // Algorithm: p1->p0 gives the vector m,
  //            n = m - (m,t) t / (t,t).
  // ----------------------------------------------------------------
  void get_n_from_t( const double &tx, const double &ty, const double &tz,
      const double &p0_x, const double &p0_y, const double &p0_z,
      const double &p1_x, const double &p1_y, const double &p1_z,
      double &nx, double &ny, double &nz );


  // ----------------------------------------------------------------
  // Calculate the circumscribing sphere's centre point and radius
  // of four given points
  // ----------------------------------------------------------------
  void get_tet_sphere_info( const double &x0, const double &x1,
      const double &x2, const double &x3, const double &y0, 
      const double &y1, const double &y2, const double &y3,
      const double &z0, const double &z1, const double &z2, 
      const double &z3, double &x, double &y, double &z, double &r );

  Vector_3 get_tet_sphere_info( const Vector_3 &pt0, const Vector_3 &pt1, 
      const Vector_3 &pt2, const Vector_3 &pt3, double &radius );

  inline double get_circumradius( const std::array<Vector_3, 4> &pts )
  {
    Matrix_double_3by3_Array AA(
        2.0 * (pts[1].x()-pts[0].x()), 2.0 * (pts[1].y()-pts[0].y()), 2.0 * (pts[1].z()-pts[0].z()),
        2.0 * (pts[2].x()-pts[0].x()), 2.0 * (pts[2].y()-pts[0].y()), 2.0 * (pts[2].z()-pts[0].z()),
        2.0 * (pts[3].x()-pts[0].x()), 2.0 * (pts[3].y()-pts[0].y()), 2.0 * (pts[3].z()-pts[0].z()) );

    AA.LU_fac();

    const double xyz2 = pts[0].dot_product( pts[0] );

    const Vector_3 centre = AA.LU_solve( Vector_3( pts[1].dot_product(pts[1]) - xyz2,
          pts[2].dot_product(pts[2]) - xyz2, pts[3].dot_product(pts[3]) - xyz2 ) );

    return ( centre - pts[0] ).norm2();
  }

  // ----------------------------------------------------------------
  // Statistical quantities
  // Mean value
  // ----------------------------------------------------------------
  inline double get_mean( const std::vector<double> &vec )
  {
    double sum = 0.0; double nn = 0.0;
    for(unsigned int ii=0; ii<vec.size(); ++ii)
    {
      sum += vec[ii];
      nn  += 1.0;
    }
    return sum / nn;
  }

  // ----------------------------------------------------------------
  // Standard deviation
  // ----------------------------------------------------------------
  inline double get_std_dev( const std::vector<double> &vec )
  {
    const double mean_val = MATH_T::get_mean(vec);
    double sum = 0.0; double nn = 0.0;
    for(unsigned int ii=0; ii<vec.size(); ++ii)
    {
      sum += (vec[ii] - mean_val) * (vec[ii] - mean_val);
      nn  += 1.0;
    }
    return std::sqrt( sum / nn );
  }
  
  // ----------------------------------------------------------------
  // Generate a Gaussian/normal distribution random value with mean value
  // mean, and standard deviation dev.
  // ----------------------------------------------------------------
  inline double gen_double_rand_normal( const double &mean, const double &std )
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::normal_distribution<double> dis(mean, std);
    return dis(gen);
  }

  // ----------------------------------------------------------------
  // gen_int_rand and gen_double_rand
  // Generate a random double in [min, max] domain for integer and double
  // precision numbers, respectively. 
  // ----------------------------------------------------------------
  inline double gen_int_rand( const int &min, const int &max )
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
  }

  inline double gen_double_rand( const double &min = -1.0, const double &max = 1.0 )
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(min, max);
    return dis(gen);
  }
  
  // ----------------------------------------------------------------
  // Print Histogram of an array of random vector
  // ----------------------------------------------------------------
  inline void print_Histogram( const std::vector<double> &val )
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

  // --------------------------------------------------------------------------
  // Projection operator
  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise constant (DGP0 space)
  // f : the function f's value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the
  //        quadrature weights.
  // nqp : number of quadrature points
  // return a scalar Prof(f) := int_omega f dx / int_omega 1 dx
  //                          = sum(f * gwts) / sum(gwts)
  // --------------------------------------------------------------------------
  double L2Proj_DGP0( const double * const &f, 
      const double * const &gwts, const int &nqp );

  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise linear (DGP1 space) in 2D.
  // f : the function value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the weights
  // qp_x : the quadrature points x-coordinates
  // qp_y : the quadrature points y-coordinates
  // nqp : the number of quadrature points
  // output: coeff_0, coeff_x, coeff_y.
  // The projected polynomial is
  //         coeff_0 + coeff_x x + coeff_y y.
  // --------------------------------------------------------------------------
  void L2Proj_DGP1_2D( const double * const &f,
      const double * const &gwts,
      const double * const &qp_x,
      const double * const &qp_y,
      const int &nqp,
      double &coeff_0, double &coeff_x, double &coeff_y );

  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise linear (DGP1 space) in 3D.
  // f : the function value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the weights
  // qp_x : the quadrature points x-coordinates
  // qp_y : the quadrature points y-coordinates
  // qp_z : the quadrature points z-coordinates
  // nqp : the number of quadrature points
  // output: coeff_0, coeff_x, coeff_y, coeff_z.
  // The projected polynomial is
  //         coeff_0 + coeff_x x + coeff_y y + coeff_z z.
  // --------------------------------------------------------------------------
  void L2Proj_DGP1_3D( const double * const &f,
      const double * const &gwts,
      const double * const &qp_x,
      const double * const &qp_y,
      const double * const &qp_z,
      const int &nqp,
      double &coeff_0, double &coeff_x, double &coeff_y, double &coeff_z );

  // ================================================================
  // Dense Matrix tool
  // This is an implementation of dense matrix in C++.
  // The matrix has to be a square matrix with size N X N.
  // The objective is to implement efficient, dense linear algebra,
  // especially the LU factorization and LU-Solver.
  // The matrix data is stored as a one-dimensional array in row-
  // oriented manner.
  // 
  // A typical usage is, one should call
  // LU_fac();
  // LU_solve(b,x);
  //
  // Ref. Numerical Linear Algebra by L.N. Trefethen and D. Bau, III, SIAM.
  // ================================================================
  class Matrix_dense
  {
    public:
      // default constructor will generate a NULL object for all ptrs
      Matrix_dense();

      // generate a in_msize x in_msize identity matrix
      Matrix_dense(const int &in_msize);

      // copy constructor
      Matrix_dense( const Matrix_dense * const & in_m );

      virtual ~Matrix_dense();

      // print the full matrix on screen
      virtual void print_mat(const int &pre = 6) const;

      // print the mat content as well as the p and invm array
      virtual void print_info() const;

      // get the size of the matrix
      virtual int get_size() const;

      virtual double get_mat(const int &ii, const int &jj) const;

      virtual double get_mat(const int &ii) const {return mat[ii];}

      virtual int get_p(const int &ii) const;

      virtual double get_invm(const int &ii) const;

      // set all entries in matrix to be zero
      virtual void zero_mat();

      // set value at i, j
      virtual void set_value(const int &ii, const int &jj, const double &val);

      virtual void set_values( int const * const &index_i, 
          int const * const &index_j,
          double const * const &vals, const int &num );

      virtual void set_values(double const * const &vals);

      // generate identity matrix
      virtual void gen_id();

      // generate matrix with random entries
      virtual void gen_rand();

      // generate Hilbert matrices
      virtual void gen_hilb();

      // multiply the mat with a vector
      // b and x canNOT be the same vector.
      virtual void Axb( double const * const &b, double * const &x ) const;

      // ------------------------------------------------------------
      // perform LU-factorization for the matrix. The mat object will be replaced
      // by the LU matrices. Only partial pivoting is performed. Complete pivoting
      // is not used because the improvement of stability is marginal and the
      // amount of time needed will increase.
      virtual void LU_fac();

      // with LU factorization performed, solve a linear problem with given RHS
      // users are responsible for allocating the b and x arrays.
      virtual void LU_solve( double const * const &b, double * const &x ) const;
      // ------------------------------------------------------------

    protected:
      // Matrix size: N by N
      const int N;

      // Length of the mat array
      const int NN;

      // double array as a holder for the matrix
      // size : N^2
      double * mat;

      // permutation infomation generated from LU-fac
      // size : N
      int * p;

      // inverse of diagonal entries which are used for LU-solve
      // size : N
      double * invm;

      // bool variable indicate if the matrix has been LU factorized.
      bool is_fac;
  };

  template<unsigned int N> class Matrix_Dense
  {
    public:
      Matrix_Dense()
      {
        for(unsigned int ii=0; ii < N*N; ++ii) mat[ii] = 0.0;
        for(unsigned int ii=0; ii < N; ++ii)
        {
          mat[ii*N + ii] = 1.0;
          p[ii] = ii;
        }
        is_fac = false;
      }

      Matrix_Dense( const std::array<double,N*N> &input )
      {
        for(unsigned int ii=0; ii < N*N; ++ii) mat[ii] = input[ii];
        for(unsigned int ii=0; ii < N; ++ii) p[ii] = ii;
        is_fac = false;
      }

      virtual ~Matrix_Dense() {};

      void print_info() const
      {
        std::cout<<"N ="<<N<<'\n';
        std::cout<<"Matrix :\n";
        int counter = -1;
        for(unsigned int ii=0; ii<N; ++ii)
        {
          for(unsigned int jj=0; jj<N; ++jj)
          {
            counter += 1;
            std::cout<<mat[counter]<<'\t';
          }
          std::cout<<'\n';
        }
        std::cout<<"p : \n";
        for(unsigned int ii=0; ii<N; ++ii)
          std::cout<<p[ii]<<'\t';
        std::cout<<std::endl;

        if(is_fac)
          std::cout<<"Matrix is factorized.\n";
        else
          std::cout<<"Matrix is NOT factorized.\n";
      }
  
      void gen_rand()
      {
        srand(time(NULL));

        for(unsigned int ii=0; ii<(N*N); ++ii)
        {
          double value = rand() % 1000; 

          mat[ii] = value * 1.0e-2 - 5;
        }

        for(unsigned int ii=0; ii<N; ++ii) p[ii] = ii;
      }

      int get_p(const int &ii) const { return p[ii];}

      double& operator()(const unsigned int &index) {return mat[index];}

      const double& operator()(const unsigned int &index) const {return mat[index];}

      double& operator()(const unsigned int &ii, const unsigned int &jj) {return mat[N*ii+jj];}

      const double& operator()(const unsigned int &ii, const unsigned int &jj) const {return mat[N*ii+jj];}

      int get_size() const {return static_cast<int>(N);}

      bool get_is_face() const {return is_fac;}

      void LU_fac()
      {
        for(unsigned int kk=0; kk<N-1; ++kk)
        {
          double max_value = std::abs(mat[kk*N+kk]);
          int max_index = kk;
          bool pivot_flag = false;
          for(unsigned int ii=kk+1; ii<N; ++ii)
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
            const int int_temp = p[kk];
            p[kk] = p[max_index];
            p[max_index] = int_temp;

            for(unsigned int ii=0; ii<N; ++ii)
            {
              const double temp = mat[kk*N+ii];
              mat[kk*N+ii] = mat[max_index*N+ii];
              mat[max_index*N+ii] = temp;
            }
          }

          const double invAkk = 1.0 / mat[kk*N+kk];

          for(unsigned int ii=kk+1; ii<N; ++ii)
          {
            mat[ii*N+kk] = mat[ii*N+kk] * invAkk;
            for(unsigned int jj=kk+1; jj<N; ++jj)
              mat[ii*N+jj] -= mat[ii*N+kk] * mat[kk*N+jj];
          }
        }

        is_fac = true;
      }

      std::array<double, N> LU_solve( std::array<double, N> &b ) const
      {
        std::array<double, N> x;
        for(unsigned int ii=0; ii<N; ++ii) x[ii] = b[p[ii]];

        for(unsigned int ii=1; ii<N; ++ii)
          for(unsigned int jj=0; jj<ii; ++jj)
            x[ii] -= mat[ii*N+jj] * x[jj];

        for(int ii=N-1; static_cast<int>(ii)>=0; --ii)
        {
          for(int jj=N-1; jj>ii; --jj)
              x[ii] -= mat[ii*N+jj] * x[jj];

          x[ii] = x[ii] / mat[ii*N+ii];
        }

        return x;
      }

      std::array<double,N> operator*( const std::array<double,N> &right )
      {
        if( is_fac == true ) std::cout<<"Warning: the matrix has been factroized.\n";
        std::array<double,N> out;
        for(unsigned int ii=0; ii<N; ++ii)
        {
          out[ii] = 0.0;
          for(unsigned int jj=0; jj<N; ++jj)
          out[ii] += mat[N*ii+jj] * right[jj];
        }
        return out;
      }

    private:
      double mat[N*N];

      unsigned int p[N];

      bool is_fac;
  };


  template<unsigned int N> std::array<double,N> operator*( const Matrix_Dense<N> &left, const std::vector<double> &right )
  {
    if( left.get_is_face() == true ) std::cout<<"Warning: the matrix has been factroized.\n";
    std::array<double,N> out;
    for(unsigned int ii=0; ii<N; ++ii)
    {
      out[ii] = 0.0;
      for(unsigned int jj=0; jj<N; ++jj)
        out[ii] += left(ii,jj) * right[jj];
    }
    return out;
  }
/*
  template<unsigned int N> std::array<double,N> operator*( const Matrix_Dense<N> &left, const std::array<double,N> &right )
  {
    if( left.get_is_face() == true ) std::cout<<"Warning: the matrix has been factroized.\n";
    std::array<double,N> out;
    for(unsigned int ii=0; ii<N; ++ii)
    {
      out[ii] = 0.0;
      for(unsigned int jj=0; jj<N; ++jj)
        out[ii] += left(ii,jj) * right[jj];
    }
    return out;
  }
*/
  template<unsigned int N> Matrix_Dense<N> operator*( const Matrix_Dense<N> &left, const Matrix_Dense<N> &right )
  {
    Matrix_Dense<N> out{};

    for(unsigned int ii=0; ii<N; ++ii)
    {
      out(ii, ii) = 0.0;

      for(unsigned int jj=0; jj<N; ++jj)
      {
        for(unsigned int kk=0; kk<N; ++kk)
        {
          out(ii, jj) += left(ii,kk) * right(kk, jj);
        }
      }
    }
    return out;
  }

  template<unsigned int N> class Matrix_SymPos_Dense : public Matrix_Dense <N>
  {
    public:
      Matrix_SymPos_Dense() : Matrix_Dense<N>()
      {}

      Matrix_SymPos_Dense(const std::array<double,N*N> &input) : Matrix_Dense<N>(&input)
      {}

      virtual ~Matrix_SymPos_Dense();

      // Check the symmetry of the matrix, throw an error if
      // non-symmetriness is found.
      void check_symm() const
      {
        for(int ii = 0; ii<N; ++ii)
        {
          for(int jj = 0; jj<ii; ++jj)
          {
            if( !MATH_T::equals( Matrix_Dense <N>::mat[ii*N+jj], Matrix_Dense <N>::mat[jj*N+ii], 1.0e-15) ) 
              std::cout<<"error: Matrix_SymPos entry ("<<ii<<","<<jj<<") does not match entry (jj"<<","<<ii<<"). \n";
          }
        }
      }
    
      // Copy the content of a matrix
      // We assume that the input matrix and the object have the same
      // size.
      void copy( const Matrix_SymPos_Dense * const &in_mat )
      {
        for(int ii=0; ii<N*N; ++ii) Matrix_Dense <N>::mat[ii] = in_mat->operator()(ii);

        for(int ii=0; ii<N; ++ii)
        {
          Matrix_Dense <N>::p[ii] = in_mat->p(ii);
        }
        Matrix_Dense <N>::is_fac = false; 
      }

      // ------------------------------------------------------------
      // Perform LDL^t transformation. The mat object will be replace
      // by the entries of the L matrix and the D matrix. Pivoting is
      // not used because this decomposition for symmetry positive
      // definite matrix is stable.
      void LDLt_fac()
      {
        // This algorithm is given in Shufang XU's book, pp 31. 
        double v[N]; 

        for(int jj=0; jj<N; ++jj)
        {
          const int Njj = jj * N;
          for(int kk=0; kk<jj; ++kk) v[kk] = Matrix_Dense <N>::mat[Njj+kk] * Matrix_Dense <N>::mat[kk*N+kk];

          for(int kk=0; kk<jj; ++kk) Matrix_Dense <N>::mat[Njj+jj] -= v[kk] * Matrix_Dense <N>::mat[Njj+kk];

          for(int ii=jj+1; ii<N; ++ii)
          {
            for(int kk=0; kk<jj; ++kk) Matrix_Dense <N>::mat[N*ii+jj] -= Matrix_Dense <N>::mat[ii*N+kk] * v[kk];

            Matrix_Dense <N>::mat[N*ii+jj] *= 1.0 / Matrix_Dense <N>::mat[Njj+jj];
          }
        }

        Matrix_Dense <N>::is_fac = true;
      }

      // With the LDLt_fac() function performed, solve a linear problem
      // with the given RHS.
      // users are responsible for allocating the b and x arrays.
      void LDLt_solve( double const * const &b, 
          double * const &x ) const
      {
        // Solve for Ly = b
        for(int ii=0; ii<N; ++ii)
        {
          x[ii] = b[ii];
          for(int jj=0; jj<ii; ++jj) x[ii] -= Matrix_Dense <N>::mat[ii*N+jj] * x[jj];
        }

        // Solve for D z = y;
        for(int ii=0; ii<N; ++ii) x[ii] *= Matrix_Dense <N>::mat[ii*N+ii];

        // Solve L^t x = z
        for(int ii=N-2; ii>=0; --ii)
        {
          for(int jj=ii+1; jj<N; ++jj) x[ii] -= Matrix_Dense <N>::mat[jj*N+ii] * x[jj];
        }
      }

      // ------------------------------------------------------------
  };


  // ================================================================
  // Dense Symmetric Positive definite matrix tool.
  // The user should be sure that the matrix is symmetric and positive
  // definite.
  // The objective is to implement efficient LDL^t decomposition,
  // which can be used to solve problems like inverting the normal
  // equation.
  // Typical usage is that one should call
  // LDLt_fac();
  // LDLt_solve(b,x);
  //
  // Ref. Shufang Xu, Numerical Linear Algebra, Peking Univ.
  // ================================================================
  class Matrix_SymPos_dense : public Matrix_dense
  {
    public:
      // Default constructor, generate a NULL object for everything
      Matrix_SymPos_dense();

      // generate a msize x msize identity matrix
      Matrix_SymPos_dense( const int &msize );

      // Copy constructor
      Matrix_SymPos_dense( const Matrix_SymPos_dense * const &in_mat );

      virtual ~Matrix_SymPos_dense();

      // Check the symmetry of the matrix, throw an error if
      // non-symmetriness is found.
      virtual void check_symm() const;

      // Copy the content of a matrix
      // We assume that the input matrix and the object have the same
      // size.
      virtual void copy( const Matrix_SymPos_dense * const &in_mat );

      // ------------------------------------------------------------
      // Perform LDL^t transformation. The mat object will be replace
      // by the entries of the L matrix and the D matrix. Pivoting is
      // not used because this decomposition for symmetry positive
      // definite matrix is stable.
      virtual void LDLt_fac();

      // With the LDLt_fac() function performed, solve a linear problem
      // with the given RHS.
      // users are responsible for allocating the b and x arrays.
      virtual void LDLt_solve( double const * const &b, 
          double * const &x ) const;
      // ------------------------------------------------------------
  };
}

#endif
