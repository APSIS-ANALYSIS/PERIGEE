#ifndef SYMMMATRIX_3X3_HPP
#define SYMMMATRIX_3X3_HPP
// ============================================================================
// SymmMatrix_3x3.hpp
// This is a 3 x 3 Symmetric Matrix class. The components are stored in a 1-D
// array following the Voigt notation. The matrix is thus
// 
//                    mat[0], mat[5], mat[4]
//                    mat[5], mat[1], mat[3]
//                    mat[4], mat[3], mat[2]
//
// Author: Yujie Sun, Ju Liu
// Date: Sept. 16 2022
// ============================================================================
#include "Matrix_3x3.hpp"

class SymmMatrix_3x3
{
  public:
    // Constructor (default an identity 3-by-3 matrix)
    SymmMatrix_3x3();

    // Copy constructor
    SymmMatrix_3x3( const SymmMatrix_3x3 &source );

    // Constructor by six numbers in Voigt numbering
    SymmMatrix_3x3( const double &m0, const double &m1, const double &m2,
        const double &m3, const double &m4, const double &m5 );

    // Destructor
    ~SymmMatrix_3x3();

    // Assignment operator
    SymmMatrix_3x3& operator= (const SymmMatrix_3x3 &source);

    // Parenthesis operator. It allows accessing and assigning the matrix entries.
    double& operator()(const int &index) {return mat[index];}

    const double& operator()(const int &index) const {return mat[index];}

    // Get-functions that access components directly
    const double& xx() const {return mat[0];}
    double& xx() {return mat[0];}

    const double& xy() const {return mat[5];}
    double& xy() {return mat[5];}

    const double& xz() const {return mat[4];}
    double& xz() {return mat[4];}

    const double& yx() const {return mat[5];}
    double& yx() {return mat[5];}

    const double& yy() const {return mat[1];}
    double& yy() {return mat[1];}

    const double& yz() const {return mat[3];}
    double& yz() {return mat[3];}

    const double& zx() const {return mat[4];}
    double& zx() {return mat[4];}

    const double& zy() const {return mat[3];}
    double& zy() {return mat[3];}

    const double& zz() const {return mat[2];}
    double& zz() {return mat[2];}
    
    // Addition operator : return left + right
    friend SymmMatrix_3x3 operator+( const SymmMatrix_3x3 &left, const SymmMatrix_3x3 &right );

    // Minus operator : return left - right
    friend SymmMatrix_3x3 operator-( const SymmMatrix_3x3 &left, const SymmMatrix_3x3 &right );

    // Add the source matrix to the object
    SymmMatrix_3x3& operator+=( const SymmMatrix_3x3 &source );

    // Minus the source matrix to the object
    SymmMatrix_3x3& operator-=( const SymmMatrix_3x3 &source );

    // Scalar product
    SymmMatrix_3x3& operator*=( const double &val );

    // Return true if the input matrix is identical to the mat
    bool is_identical( const SymmMatrix_3x3 &source, const double &tol = 1.0e-12 ) const;

    // Set all components to zero
    void gen_zero();

    // Set an identity matrix
    void gen_id();

    // Set components a random value
    void gen_rand();

    // Invert the matrix
    void inverse();

    // add the matrix with a given matrix with scaling
    // X = X + a * Y
    void AXPY( const double &val, const SymmMatrix_3x3 &source );

    // X = X + a * I
    void AXPI( const double &val );

    // Get the determinant of the matrix
    double det() const;

    // Get the trace of the matrix
    double tr() const {return mat[0] + mat[1] + mat[2];}

    // Get the invariants
    double I1() const {return tr();}

    double I2() const;

    double I3() const {return det();}

    // Return x^T Mat y, assuming x, y are both column vectors of size 3
    double VecMatVec( const Vector_3 &x, const Vector_3 &y ) const;

    // Vector multiplication y = Ax, the vectors have to be size 3
    Vector_3 VecMult( const Vector_3 &x ) const;

    void VecMult( const double &x0, const double &x1, const double &x2,
       double &y0, double &y1, double &y2 ) const;

    // Matrix rotation
    // Let Q be a rotation matrix, the matrix gets updated by
    // Q^T M Q = Q_ki M_kl Q_lj = output_matrix_ij
    void MatRot( const Matrix_3x3 &Q );

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
    
    double MatContraction( const SymmMatrix_3x3 &source ) const;

    // print the matrix
    void print() const;

    // print the matrix in a row
    void print_in_row() const;

    // print the Voigt components in the order of xx yy zz yz xz xy
    void print_Voigt() const;

  private:
    double mat[6];
};

Vector_3 operator*( const SymmMatrix_3x3 &left, const Vector_3 &right );

Matrix_3x3 operator*( const SymmMatrix_3x3 &left, const Matrix_3x3 &right );

Matrix_3x3 operator*( const Matrix_3x3 &left, const SymmMatrix_3x3 &right );

Matrix_3x3 operator*( const SymmMatrix_3x3 &left, const SymmMatrix_3x3 &right );

SymmMatrix_3x3 operator*( const double &val, const SymmMatrix_3x3 &input );

// Return the inverse of the input matrix
SymmMatrix_3x3 inverse( const SymmMatrix_3x3 &input );

// C = F^T F
SymmMatrix_3x3 gen_right_Cauchy_Green( const Matrix_3x3 &input );

// b = F F^T
SymmMatrix_3x3 gen_left_Cauchy_Green( const Matrix_3x3 &input );

// convert a regular matrix to its symmetric part
// output = 0.5 x ( source + source_transpose )
SymmMatrix_3x3 gen_symm_part( const Matrix_3x3 &input );

#endif
