#ifndef MATRIX_DENSE_HPP
#define MATRIX_DENSE_HPP
// ============================================================================
// Matrix_dense.hpp
// This is an implementation of dense matrix in C++.
// 
// This matrix has to be a square matrix. Its size is N x N.
//
// The objective is to implement fast dense linear algebra for NxN dense
// matrices. In particular, LU factorization and LU-solver are implemented.
//
// The date are saved in a one-dimensional array for the matrix in row-oriented
// manner.
//
// 
// Reference: Numerical Linear Algebra by L.N. Trefethen and D. Bau, III, SIAM.
//
// Date: Oct. 1st 2015
// ============================================================================
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

class Matrix_dense
{
  public:
    // default constructor will generate a NULL object for all ptrs
    Matrix_dense();

    // generate a in_msize x in_msize identity matrix
    Matrix_dense(const int &in_msize);

    // copy constructor
    Matrix_dense( const Matrix_dense * const & in_m );

    ~Matrix_dense();

    // print the full matrix on screen
    void print_mat(const int &pre = 6) const;

    // print the mat content as well as the p and invm array
    void print_info() const;

    // get the size of the matrix
    int get_size() const;

    double get_mat(const int &ii, const int &jj) const;
    
    int get_p(const int &ii) const;

    double get_invm(const int &ii) const;

    // set all entries in matrix to be zero
    void zero_mat();

    // set value at i, j
    void set_value(const int &ii, const int &jj, const double &val);

    void set_values( int const * const &index_i, int const * const &index_j,
        double const * const &vals, const int &num );

    void set_values(double const * const &vals);

    // generate matrix with random entries
    void gen_rand();
    
    // generate Hilbert matrices
    void gen_hilb();

    // perform LU-factorization for the matrix. The mat object will be replaced
    // by the LU matrices. Only partial pivoting is performed. Complete pivoting
    // is not used because the improvement of stability is marginal and the
    // amout of time needed will increase.
    void LU_fac();
    
    // with LU factorization performed, solve a linear problem with givne RHS
    void LU_solve( double const * const &b, double * const &x ) const;

    // multiply the mat with a vector
    void Axb( double const * const &b, double * const &x ) const;

  private:
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
    bool is_LU;
};

#endif
