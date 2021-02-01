#ifndef MATRIX_3X3_HPP
#define MATRIX_3X3_HPP
// ==================================================================
// Matrix_3x3.hpp
// This is a 3 x 3 matrix class. The components are stored in a 1-D 
// array: mat[9]. Logically, the matrix is 
//
//                   mat[0], mat[1], mat[2]
//                   mat[3], mat[4], mat[5]
//                   mat[6], mat[7], mat[8]
//
// Note: This is adopted from a previous implementation Matrix_double
// _3by3_Array.hpp. The difference is that we do not intend to 
// implement the LU factorization in this class; this means that we do
// not have the pivoting flag array and the inverse diagonal array in 
// this class. If needed, one can implement a LU-factorization in this
// class with pivoting and LU factorized matrix as an output.
//
// Author: Ju Liu
// Date: June 21 2016
// ==================================================================
#include "Vector_3.hpp"

class Matrix_3x3
{
  public:
    // Constructor (default an identity 3-by-3 matrix)
    Matrix_3x3();

    // Copy constructor
    Matrix_3x3( const Matrix_3x3 &source );

    // Explicit Defintion of all 9 entries
    Matrix_3x3( const double &a11, const double &a12, const double &a13,
        const double &a21, const double &a22, const double &a23,
        const double &a31, const double &a32, const double &a33 );

    // Destructor
    ~Matrix_3x3();

    // Copy
    void copy( const Matrix_3x3 &source );    
    
    void copy( double source[9] );

    // Assignment operator
    Matrix_3x3& operator= (const Matrix_3x3 &source);

    // Parenthesis operator. It allows accessing and assigning the matrix
    // entries.
    double& operator()(const int &index) {return mat[index];}

    const double& operator()(const int &index) const {return mat[index];}

    // Parenthesis operator. Access through row and col index: ii jj
    // Note: We do not check that ii , jj = 0, 1, 2.
    double& operator()(const int &ii, const int &jj)
    {return mat[3*ii+jj];}

    const double& operator()(const int &ii, const int &jj) const
    {return mat[3*ii+jj];}

    // Addition operator : return left + right
    friend Matrix_3x3 operator+( const Matrix_3x3 &left, const Matrix_3x3 &right);

    // Minus operator : return left - right
    friend Matrix_3x3 operator-( const Matrix_3x3 &left, const Matrix_3x3 &right);
    
    // Add the source matrix to the object
    Matrix_3x3& operator+=( const Matrix_3x3 &source );

    // Minus the source matrix to the object
    Matrix_3x3& operator-=( const Matrix_3x3 &source );

    // Scalar product
    Matrix_3x3& operator*=( const double &val );

    // Return true if the input matrix is identical to the mat
    bool is_identical( const Matrix_3x3 source ) const;

    // Set all components to zero
    void gen_zero();

    // Set an identity matrix
    void gen_id();

    // Set components a random value
    void gen_rand();

    // Set a Hilbert matrix
    void gen_hilb();

    // Set a matrix from out-product of two vecs with length 3
    // mat_ij = a_i b_j
    void gen_outprod( const double * const &a, const double * const &b );
    
    void gen_outprod( const Vector_3 &a, const Vector_3 &b );
   
    // mat_ij = a_i a_j 
    void gen_outprod( const double * const &a );
    
    void gen_outprod( const Vector_3 &a );

    // Transpose the matrix
    void transpose();

    // Invert the matrix
    void inverse();

    // Scale the matrix by a scalar
    void scale( const double &val );

    // add the matrix with a given matrix with scaling
    // X = X + a * Y
    void AXPY( const double &val, const Matrix_3x3 &source );

    // add the matrix source with the matrix
    // X = X + Y
    void PY( const Matrix_3x3 &source );

    // Get the determinant of the matrix
    double det() const;

    // Get the trace of the matrix
    double tr() const {return mat[0] + mat[4] + mat[8];}

    // Get the invariants
    double I1() const {return tr();}

    double I2() const;

    double I3() const {return det();}

    // Return x^T Mat y, assuming x, y are both column vectors of size 3
    double VecMatVec( const double * const &x, const double * const &y ) const;

    double VecMatVec( const Vector_3 &x, const Vector_3 &y ) const;

    // Vector multiplication y = Ax, the vectors have to be size 3
    void VecMult( const double * const &x, double * const &y ) const;
    
    void VecMult( const Vector_3 &x, Vector_3 &y ) const;

    // y = Ax, wherein x = [x0; x1; x2]
    void VecMult( const double &x0, const double &x1, const double &x2,
       double * const &y ) const;

    // y = x^T A,  in indices: y_i = x_I A_Ii
    void VecMultT( const double * const &x, double * const &y ) const;
    
    void VecMultT( const Vector_3 &x, Vector_3 &y ) const;

    // y = x^T A, wherein x = [x0; x1; x2]
    void VecMultT( const double &x0, const double &x1, const double &x2,
       double * const &y ) const;

    // x = Ax, vector x has to be of size 3
    void VecMult( double * const &x ) const;
    
    void VecMult( Vector_3 &x ) const;

    // Matrix multiplication mat = mleft * mright
    void MatMult( const Matrix_3x3 &mleft, const Matrix_3x3 &mright );

    // Matrix multiplication as mat = source^T * source
    // This is used for the evaluation of right Cauchy-Green strain tensor:
    //                       C = F^T F
    // The resulting matrix is symmetric. Hence the computation is simplified.
    void MatMultTransposeLeft( const Matrix_3x3 &source );

    // Matrix multiplication as mat = source * source^T
    // This is used for the evaluation of the left Cauchy-Green strain tensor:
    //                       b = F F^T
    // The resulting matrix is symmetric. Hence, the computation is simplified.
    void MatMultTransposeRight( const Matrix_3x3 &source );

    // Matrix contraction
    // return mat_ij source_ij
    double MatContraction( const Matrix_3x3 &source ) const;

    // return mat_ij source_ji
    double MatTContraction( const Matrix_3x3 &source ) const;

    // print the matrix
    void print() const;

    // Eigen decomposition of the matrix M = eta1 v1 v1T + eta2 v2 v2T + eta3 v3
    // v3T. The algorithm is based on CMAME 197 2008 4007-4015 paper by
    // W.M. Scherzinger and C.R. Dohrmann
    // return 1 if the three eigenvalues are the same
    // return 2 if there are two identical eigenvalues, the most distinct one is
    // eta_1 and eta_1's associated eigenvector is v1
    // return 3 if all three are distinct
    int eigen_decomp( double &eta1, double &eta2, double &eta3,
       Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const;

  private:
    double mat[9];

    // Find the eignevector correspond to a eigenvalue
    // This implements the algorithm documented in
    // CMAME 2008 v197 4007-4015, sec 2.4
    // It will return an eigenvector v for this eigenvalue,
    // and two additional vectors s1 s2; v-s1-s2 forms a orthonormal basis.
    void find_eigen_vector( const double &eta, Vector_3 &v,
        Vector_3 &s1, Vector_3 &s2 ) const;

    // Return the deviatoric component's contraction scaled by 0.5.
    // M' = M - 0.3333 tr(M) I, return 0.5 M : M.
    double J2() const;

    // Return the determinant of the deviatoric component.
    // M' = M - 0.3333 tr(M) I, return det(M').
    double J3() const;
};

#endif
