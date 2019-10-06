#ifndef MATRIX_2X2_HPP
#define MATRIX_2X2_HPP
// ==================================================================
// Matrix_2x2.hpp
// This is a 2 x 2 matrix class. The components are stored in a 1-D
// array: mat[4], or logically
//                  mat[0], mat[1]
//                  mat[2], mat[3]
// The class is designed for the 2nd-order tensor in 2D.
//
// Author: Ju Liu
// Date: Sept. 12 2016
// ==================================================================
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

class Matrix_2x2
{
  public:
    // Default Constructor: generate an identity 2-by-2 matrix
    Matrix_2x2();

    // Copy constructor
    Matrix_2x2( const Matrix_2x2 &source );

    // Explicit defintion of all 4 components
    Matrix_2x2( const double &a11, const double &a12,
        const double &a21, const double &a22 );

    // Destructor
    ~Matrix_2x2();

    // Copy
    void copy( const Matrix_2x2 &source );

    void copy( double source[4] );

    // Parenthesis operators.
    double& operator()(const int &index) {return mat[index];}

    const double& operator()(const int &index) const {return mat[index];}

    double& operator()(const int &ii, const int &jj)
    {return mat[2*ii+jj];}

    const double& operator()(const int &ii, const int &jj) const
    {return mat[2*ii+jj];}

    // Set all components to zero
    void gen_zero();

    // Set an identity matrix
    void gen_id();

    // Set components a random value
    void gen_rand();

    // Set a Hilbert matrix
    void gen_hilb();

    // Transpose the matrix
    void transpose();

    // Invert the matrix
    void inverse();

    // Scale the matrix by a scalar
    void scale( const double &val );

    // add the matrix with a given matrix with scaling
    // X = X + a * Y
    void AXPY( const double &val, const Matrix_2x2 &source );

    // add the matrix source with the matrix
    // X = X + Y
    void PY( const Matrix_2x2 &source );

    // Get the determinant of the matrix
    double det() const;

    // Get the trace of the matrix
    double tr() const {return mat[0] + mat[3];}

    // Return x^T Mat y, assuming x, y are both column vectors of size 2
    double VecMatVec( const double * const &x, 
        const double * const &y ) const;

    // Vector multiplication y = Ax, the vectors have to be size 2
    void VecMult( const double * const &x, double * const &y ) const;

    // x = Ax, vector x has to be of size 2
    void VecMult( double * const &x ) const;

    // Matrix multiplication mat = mleft * mright
    void MatMult( const Matrix_2x2 &mleft, const Matrix_2x2 &mright );


    // Matrix multiplication as mat = source^T * source
    // This is used for the evaluation of right Cauchy-Green strain tensor:
    //                       C = F^T F
    // The resulting matrix is symmetric. Hence the computation is simplified.
    void MatMultTransposeLeft( const Matrix_2x2 &source );

    // Matrix multiplication as mat = source * source^T
    // This is used for the evaluation of the left Cauchy-Green strain tensor:
    //                       b = F F^T
    // The resulting matrix is symmetric. Hence, the computation is simplified.
    void MatMultTransposeRight( const Matrix_2x2 &source );

    // Matrix contraction
    // return mat_ij source_ij
    double MatContraction( const Matrix_2x2 &source ) const;

    // print the matrix
    void print() const;

  private:
    double mat[4];
};

#endif
