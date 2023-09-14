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
#include <random>
#include <iostream>
#include <iomanip>

namespace MATH_T
{
  // --------------------------------------------------------------------------
  // Useful Constants:
  // --------------------------------------------------------------------------
  // PI = 3.1415926
  const double PI = atan(1.0) * 4.0;

  // E = 2.71828
  const double E  = exp(1.0);

  // Compute the binomial coefficient, e.g.
  // n = 0 : 1                  n = 1 : 1 1
  // n = 2 : 1 2 1              n = 3 : 1 3 3 1
  // n = 4 : 1 4 6 4 1          n = 5 : 1 5 10 10 5 1
  inline double binomialCoefficient( const double &n, const double &k )
  {
    if( (k<0) || (k>n) ) return 0;

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

  // --------------------------------------------------------------------------
  // This function is used to determine if two vector object are identical
  // up to a tolerance (default 1.0e-12).
  // --------------------------------------------------------------------------
  template<typename T> bool equals( const std::vector<T> &a, 
      const std::vector<T> &b, const double &tol = 1.0e-12 )
  {
    if( a.size() != b.size() ) return false;
    for(unsigned int ii=0; ii<a.size(); ++ii)
    {
      if( std::abs(a[ii]-b[ii]) >= tol ) return false;
    }
    return true;
  }
  
  // --------------------------------------------------------------------------
  // Cross product of two 3D vectors
  // Input: u = (u1, u2, u3)
  //        v = (v1, v2, v3)
  // Output: uxv = (u2v3 - u3v2)i + (u3v1 - u1v3)j + (u1v2 - u2v1)k
  // --------------------------------------------------------------------------
  inline void cross3d(const double &u1, const double &u2, const double &u3,
      const double &v1, const double &v2, const double &v3,
      double &w1, double &w2, double &w3 )
  {
    w1 = u2 * v3 - u3 * v2;
    w2 = u3 * v1 - u1 * v3;
    w3 = u1 * v2 - u2 * v1;
  }

  // --------------------------------------------------------------------------
  // Normalize 3D vector
  // Input: x, y, z
  // Output: x/len, y/len, z/len, len = sqrt(x^2+y^2+z^2)
  // --------------------------------------------------------------------------
  inline double normalize3d( double &x, double &y, double &z )
  {
    const double len = std::sqrt(x*x + y*y + z*z);
    x = x / len;
    y = y / len;
    z = z / len;

    return len;
  }

  inline double norm2(const double &x, const double &y, const double &z)
  {
    return std::sqrt(x*x + y*y + z*z);
  }

  // ----------------------------------------------------------------
  // Statistical quantities
  // Mean value
  // ----------------------------------------------------------------
  inline double get_mean( const std::vector<double> &vec )
  {
    double sum = 0.0; double nn = 0.0;
    const unsigned int len = vec.size();
    for(unsigned int ii=0; ii<len; ++ii)
    {
      sum += vec[ii];
      nn  += 1.0;
    }
    return sum / nn;
  }

  // ----------------------------------------------------------------
  // Standard deviation
  // ----------------------------------------------------------------
  double get_std_dev( const std::vector<double> &vec );
  
  // ----------------------------------------------------------------
  // Generate a Gaussian distribution vector with length n, mean value
  // mean, and standard deviation dev, using Marsaglia algorithm
  //
  // Note: Call srand((unsigned int)time(NULL)) before calling this generator!
  // ----------------------------------------------------------------
  void gen_Gaussian( const int &n, const double &mean, const double &std,
      std::vector<double> &val );

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
  void print_Histogram( const std::vector<double> &val );
  
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
