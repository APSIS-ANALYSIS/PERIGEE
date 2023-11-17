#ifndef TENSOR2_2D_HPP
#define TENSOR2_2D_HPP
// ==================================================================
// Tensor2_2D.hpp
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
#include <random>

class Tensor2_2D
{
  public:
    // Default Constructor: generate an identity 2-by-2 matrix
    Tensor2_2D();

    // Copy constructor
    Tensor2_2D( const Tensor2_2D &source );

    // Explicit defintion of all 4 components
    Tensor2_2D( const double &a11, const double &a12,
        const double &a21, const double &a22 );

    // Destructor
    ~Tensor2_2D() = default;

    // Copy
    void copy( const Tensor2_2D &source );

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
    void gen_rand(const double &min = -1.0, const double &max = 1.0);

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
    void AXPY( const double &val, const Tensor2_2D &source );

    // add the matrix source with the matrix
    // X = X + Y
    void PY( const Tensor2_2D &source );

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
    void MatMult( const Tensor2_2D &mleft, const Tensor2_2D &mright );


    // Matrix multiplication as mat = source^T * source
    // This is used for the evaluation of right Cauchy-Green strain tensor:
    //                       C = F^T F
    // The resulting matrix is symmetric. Hence the computation is simplified.
    void MatMultTransposeLeft( const Tensor2_2D &source );

    // Matrix multiplication as mat = source * source^T
    // This is used for the evaluation of the left Cauchy-Green strain tensor:
    //                       b = F F^T
    // The resulting matrix is symmetric. Hence, the computation is simplified.
    void MatMultTransposeRight( const Tensor2_2D &source );

    // Matrix contraction
    // return mat_ij source_ij
    double MatContraction( const Tensor2_2D &source ) const;

    // print the matrix
    void print() const;

  private:
    double mat[4];
};

#endif
