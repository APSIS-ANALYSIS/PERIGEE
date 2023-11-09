#ifndef TENSOR2_3D_HPP
#define TENSOR2_3D_HPP
// ============================================================================
// Tensor2_3D.hpp
// This is a 2nd-order tensor in 3D. The components are stored in a 1-D array: 
// mat[9]. Logically, the matrix is 
//
//                   mat[0], mat[1], mat[2]
//                   mat[3], mat[4], mat[5]
//                   mat[6], mat[7], mat[8]
//
// Note: We do not intend to implement the LU factorization in this class; this
// means that we do not have the pivoting flag array in this class. If one wants
// to solve a small dense matrix problem efficiently, one should refer to the
// matrix class implemented in Math_Tools.hpp. This class is used primarily for
// material modeling.
//
// Author: Ju Liu
// Date: June 21 2016
// ============================================================================
#include "Vector_3.hpp"

class Tensor2_3D
{
  public:
    // Constructor (default an identity 3-by-3 matrix)
    Tensor2_3D();

    // Copy constructor
    Tensor2_3D( const Tensor2_3D &source );

    // Explicit Defintion of all 9 entries
    Tensor2_3D( const double &a11, const double &a12, const double &a13,
        const double &a21, const double &a22, const double &a23,
        const double &a31, const double &a32, const double &a33 );

    // Generate a matrix made by 3 column vectors [ vec1 | vec2 | vec3 ]
    Tensor2_3D ( const Vector_3 &vec1, const Vector_3 &vec2, const Vector_3 &vec3 );
    
    // Destructor
    ~Tensor2_3D();

    // Copy
    void copy( const Tensor2_3D &source );    
    
    // Assignment operator
    Tensor2_3D& operator= (const Tensor2_3D &source);

    // Parenthesis operator. It allows accessing and assigning the matrix entries.
    double& operator()(const int &index) {return mat[index];}

    const double& operator()(const int &index) const {return mat[index];}

    // Parenthesis operator. Access through row and col index: ii jj
    // Note: We do not check that ii , jj = 0, 1, 2.
    double& operator()(const int &ii, const int &jj) {return mat[3*ii+jj];}

    const double& operator()(const int &ii, const int &jj) const
    {return mat[3*ii+jj];}

    // Get-functions that access components directly via the get-function's name
    const double& xx() const {return mat[0];}
    double& xx() {return mat[0];}
    
    const double& xy() const {return mat[1];}
    double& xy() {return mat[1];}
    
    const double& xz() const {return mat[2];}
    double& xz() {return mat[2];}

    const double& yx() const {return mat[3];}
    double& yx() {return mat[3];}
    
    const double& yy() const {return mat[4];}
    double& yy() {return mat[4];}
    
    const double& yz() const {return mat[5];}
    double& yz() {return mat[5];}

    const double& zx() const {return mat[6];}
    double& zx() {return mat[6];}
    
    const double& zy() const {return mat[7];}
    double& zy() {return mat[7];}
    
    const double& zz() const {return mat[8];}
    double& zz() {return mat[8];}

    // Addition operator : return left + right
    friend Tensor2_3D operator+( const Tensor2_3D &left, const Tensor2_3D &right);

    // Minus operator : return left - right
    friend Tensor2_3D operator-( const Tensor2_3D &left, const Tensor2_3D &right);
    
    // Add the source matrix to the object
    Tensor2_3D& operator+=( const Tensor2_3D &source );

    // Minus the source matrix to the object
    Tensor2_3D& operator-=( const Tensor2_3D &source );

    // Scalar product
    Tensor2_3D& operator*=( const double &val );

    // unary minus operator
    Tensor2_3D operator-() const;

    // Return true if the input matrix is identical to the mat
    bool is_identical( const Tensor2_3D &source, const double &tol = 1.0e-12 ) const;

    // Set all components to zero
    void gen_zero();

    // Set an identity matrix
    void gen_id();

    // Set components a random value
    void gen_rand(const double &left = -1.0, const double &right = 1.0);

    // Set a Hilbert matrix
    void gen_hilb();

    // Set a matrix from out-product of two vecs with length 3
    // mat_ij = a_i b_j
    void gen_outprod( const Vector_3 &va, const Vector_3 &vb );
   
    // mat_ij = a_i a_j 
    void gen_outprod( const Vector_3 &va );
    
    // Add a matrix from out-product of two vecs with length 3
    // mat_ij += val a_i b_j
    void add_outprod( const double &val, const Vector_3 &va, const Vector_3 &vb );
   
    // Transpose the matrix
    void transpose();

    // Invert the matrix
    void inverse();

    // Scale the matrix by a scalar
    void scale( const double &val );

    // add the matrix with a given matrix with scaling
    // X = X + a * Y
    void AXPY( const double &val, const Tensor2_3D &source );

    // X = X + a * I
    void AXPI( const double &val );

    // Get the determinant of the matrix
    double det() const;

    // Get the trace of the matrix
    double tr() const {return mat[0] + mat[4] + mat[8];}

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

    // y = x^T A,  in indices: y_i = x_I A_Ii
    Vector_3 VecMultT( const Vector_3 &x ) const;

    void VecMultT( const double &x0, const double &x1, const double &x2,
       double &y0, double &y1, double &y2 ) const;

    // Matrix multiplication mat = mleft * mright
    void MatMult( const Tensor2_3D &mleft, const Tensor2_3D &mright );

    // ------------------------------------------------------------------------
    // Matrix rotation
    // Let Q be a rotation matrix, the matrix gets updated by
    // Q^T M Q = Q_ki M_kl Q_lj = output_matrix_ij
    // ------------------------------------------------------------------------
    void MatRot( const Tensor2_3D &Q );

    // ------------------------------------------------------------------------
    // Matrix multiplication as mat = source^T * source
    // This is used for the evaluation of right Cauchy-Green strain tensor:
    //                       C = F^T F
    // The resulting matrix is symmetric. Hence the computation is simplified.
    // ------------------------------------------------------------------------
    void MatMultTransposeLeft( const Tensor2_3D &source );

    // ------------------------------------------------------------------------
    // Matrix multiplication as mat = source * source^T
    // This is used for the evaluation of the left Cauchy-Green strain tensor:
    //                       b = F F^T
    // The resulting matrix is symmetric. Hence, the computation is simplified.
    // ------------------------------------------------------------------------
    void MatMultTransposeRight( const Tensor2_3D &source );

    // ------------------------------------------------------------------------
    // Matrix contraction
    // return mat_ij source_ij
    // ------------------------------------------------------------------------
    double MatContraction( const Tensor2_3D &source ) const;

    // ------------------------------------------------------------------------
    // Contraction with the transposed input
    // return mat_ij source_ji
    // ------------------------------------------------------------------------
    double MatTContraction( const Tensor2_3D &source ) const;

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
    // Note: this function is applicable to symmetric tensors only. Yet, we do
    // not perform symmetry check within this function.
    // ------------------------------------------------------------------------
    int eigen_decomp( double &eta1, double &eta2, double &eta3,
       Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const;

  private:
    double mat[9];

    // ------------------------------------------------------------------------
    // Find the eignevector correspond to a eigenvalue
    // This implements the algorithm documented in
    // CMAME 2008 v197 4007-4015, sec 2.4
    // It will return an eigenvector v for this eigenvalue,
    // and two additional vectors s1 s2; v-s1-s2 forms orthonormal bases.
    // ------------------------------------------------------------------------
    void find_eigen_vector( const double &eta, Vector_3 &v,
        Vector_3 &s1, Vector_3 &s2 ) const;

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

// Return the matrix-vectror multiplication
Vector_3 operator*( const Tensor2_3D &left, const Vector_3 &right );

// Return the matrix-matrix multiplication
Tensor2_3D operator*( const Tensor2_3D &left, const Tensor2_3D &right );

// Return the scalar scaling of matrix
Tensor2_3D operator*( const double &val, const Tensor2_3D &input );

namespace Ten2
{
  // Return the inverse of the input matrix
  Tensor2_3D inverse( const Tensor2_3D &input );

  // Return the cofactor of input matrix which is J input^(-T)
  Tensor2_3D cofactor( const Tensor2_3D &input );

  // Return the transpose of input matrix
  Tensor2_3D transpose( const Tensor2_3D &input );

  // Return an identity matrix
  Tensor2_3D gen_id();

  // Return a zero matrix
  Tensor2_3D gen_zero();

  // Return the exponential of the input matrix
  // exp(X) = sum_{k=0}^{infty} 1/(k!) X^k
  // The algorithm can be found on page 749 of Computational Methods for Plasticity
  // written by E. A. De Souza Neto, D. Perić and D. R. J. Owen.
  Tensor2_3D exp( const Tensor2_3D &input );
}

#endif
