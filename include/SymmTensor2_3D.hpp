#ifndef SYMMTENSOR2_3D_HPP
#define SYMMTENSOR2_3D_HPP
// ============================================================================
// SymmTensor2_3D.hpp
// This is a 3 x 3 Symmetric Matrix/Tensor class. The components are stored in
// a 1-D array following the Voigt notation. The matrix is thus
// 
//                    mat[0], mat[5], mat[4]
//                    mat[5], mat[1], mat[3]
//                    mat[4], mat[3], mat[2]
//
// Author: Yujie Sun, Ju Liu
// Date: Sept. 16 2022
// ============================================================================
#include "Tensor2_3D.hpp"

class SymmTensor2_3D
{
  public:
    // Constructor (default an identity 3-by-3 matrix)
    SymmTensor2_3D();

    // Copy constructor
    SymmTensor2_3D( const SymmTensor2_3D &source );

    // Constructor by six numbers in Voigt numbering
    SymmTensor2_3D( const double &m0, const double &m1, const double &m2,
        const double &m3, const double &m4, const double &m5 );

    // Destructor
    ~SymmTensor2_3D() = default;

    // Assignment operator
    SymmTensor2_3D& operator= (const SymmTensor2_3D &source);

    Tensor2_3D convert_to_full() const;

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
    friend SymmTensor2_3D operator+( const SymmTensor2_3D &left, const SymmTensor2_3D &right );

    // Minus operator : return left - right
    friend SymmTensor2_3D operator-( const SymmTensor2_3D &left, const SymmTensor2_3D &right );

    // Add the source matrix to the object
    SymmTensor2_3D& operator+=( const SymmTensor2_3D &source );

    // Minus the source matrix to the object
    SymmTensor2_3D& operator-=( const SymmTensor2_3D &source );

    // Scalar product
    SymmTensor2_3D& operator*=( const double &val );

    // Return true if the input matrix is identical to the mat
    bool is_identical( const SymmTensor2_3D &source, const double &tol = 1.0e-12 ) const;

    // Set all components to zero
    void gen_zero();

    // Set an identity matrix
    void gen_id();

    // Set components a random value
    void gen_rand(const double &left = -1.0, const double &right = 1.0);

    // Invert the matrix
    void inverse();

    // add the matrix with a given matrix with scaling
    // X = X + a * Y
    void AXPY( const double &val, const SymmTensor2_3D &source );

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
    void MatRot( const Tensor2_3D &Q );

    // Matrix contraction
    // return mat_ij source_ij
    double MatContraction( const Tensor2_3D &source ) const;
    
    double MatContraction( const SymmTensor2_3D &source ) const;

    // Obtain the Voigt notation for a regular matrix index
    // 0 <= index < 9
    int Voigt_notation( const int &index ) const
    {
      // This map is used to transform the natural indices of a 3x3 symmetric matrix
      // to Voigt notation
      constexpr int map[9] = { 0, 5, 4,
        5, 1, 3,
        4, 3, 2 };
      return map[index];
    }

    // print the matrix
    void print() const;

    // print the matrix in a row
    void print_in_row() const;

    // print the Voigt components in the order of xx yy zz yz xz xy
    void print_Voigt() const;

    // ------------------------------------------------------------------------
    // Eigen decomposition of the matrix 
    // M = eta1 v1 v1T + eta2 v2 v2T + eta3 v3 v3T. 
    // The algorithm is based on CMAME 197 2008 4007-4015 paper by W.M. Scherzinger
    // and C.R. Dohrmann.
    // This function will
    // return 1 if the three eigenvalues are the same
    // return 2 if there are two identical eigenvalues, the most distinct one is
    // eta_1 and eta_1's associated eigenvector is v1
    // return 3 if all three are distinct
    // ------------------------------------------------------------------------
    int eigen_decomp( double &eta1, double &eta2, double &eta3,
       Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const;

    // ------------------------------------------------------------------------
    // Find the eignevector correspond to a eigenvalue
    // This implements the algorithm documented in
    // CMAME 2008 v197 4007-4015, sec 2.4
    // It will return an eigenvector v for this eigenvalue,
    // and two additional vectors s1 s2; v-s1-s2 forms a orthonormal basis.
    // ------------------------------------------------------------------------
    void find_eigen_vector( const double &eta, Vector_3 &v,
        Vector_3 &s1, Vector_3 &s2 ) const;

  private:
    double mat[6];

    // ------------------------------------------------------------------------
    // Return the deviatoric component's contraction scaled by 0.5.
    // M' = M - 0.3333 tr(M) I, return 0.5 M' : M'.
    // ------------------------------------------------------------------------
    double J2() const;

    // ------------------------------------------------------------------------
    // Return the determinant of the deviatoric component.
    // M' = M - 0.3333 tr(M) I, return det(M').
    // ------------------------------------------------------------------------
    double J3() const;
};

Vector_3 operator*( const SymmTensor2_3D &left, const Vector_3 &right );

Tensor2_3D operator*( const SymmTensor2_3D &left, const Tensor2_3D &right );

Tensor2_3D operator*( const Tensor2_3D &left, const SymmTensor2_3D &right );

Tensor2_3D operator*( const SymmTensor2_3D &left, const SymmTensor2_3D &right );

SymmTensor2_3D operator*( const double &val, const SymmTensor2_3D &input );

namespace STen2
{
  // Return the inverse of the input matrix
  SymmTensor2_3D inverse( const SymmTensor2_3D &input );

  // Generate the right Cauchy-Green tensor C = F^T F
  SymmTensor2_3D gen_right_Cauchy_Green( const Tensor2_3D &input );

  // Generate the left Cauchy-Green tensor b = F F^T
  SymmTensor2_3D gen_left_Cauchy_Green( const Tensor2_3D &input );

  // Convert a regular matrix to its symmetric part
  // output = 0.5 x ( source + source_transpose )
  SymmTensor2_3D gen_symm_part( const Tensor2_3D &input );
}

#endif
