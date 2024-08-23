#ifndef TENSOR2_3D_HPP
#define TENSOR2_3D_HPP
// ============================================================================
// Tensor2_3D.hpp
// This is a 2nd-order tensor in 3D. The components are stored in a 1-D array: 
// mat[9]. Logically, the tensor components are organized as 
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
    constexpr Tensor2_3D() : mat{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}} {}

    // Copy constructor
    Tensor2_3D( const Tensor2_3D &source ) : mat(source.mat) {}

    // Explicit Defintion of all 9 entries
    constexpr Tensor2_3D( const double &a11, const double &a12, const double &a13,
        const double &a21, const double &a22, const double &a23,
        const double &a31, const double &a32, const double &a33 )
      : mat{{a11, a12, a13, a21, a22, a23, a31, a32, a33}} {}

    // Generate a matrix made by 3 column vectors [ vec1 | vec2 | vec3 ]
    constexpr Tensor2_3D ( const Vector_3 &vec1, const Vector_3 &vec2, const Vector_3 &vec3 )
      : mat{{ vec1(0), vec2(0), vec3(0), vec1(1), vec2(1), vec3(1),
        vec1(2), vec2(2), vec3(2) }} {}
    
    // Destructor
    ~Tensor2_3D() = default;

    // Assignment operator
    Tensor2_3D& operator= (const Tensor2_3D &source);

    // Parenthesis operator. It allows accessing and assigning the matrix entries.
    inline double& operator()(const int &index) {return mat[index];}

    inline const double& operator()(const int &index) const {return mat[index];}

    // Parenthesis operator. Access through row and col index: ii jj
    // Note: We do not check that ii , jj = 0, 1, 2.
    inline double& operator()(const int &ii, const int &jj) 
    {return mat[3*ii+jj];}

    inline const double& operator()(const int &ii, const int &jj) const
    {return mat[3*ii+jj];}

    // Get-functions that access components directly via the get-function's name
    inline const double& xx() const {return mat[0];}
    inline double& xx() {return mat[0];}
    
    inline const double& xy() const {return mat[1];}
    inline double& xy() {return mat[1];}
    
    inline const double& xz() const {return mat[2];}
    inline double& xz() {return mat[2];}

    inline const double& yx() const {return mat[3];}
    inline double& yx() {return mat[3];}
    
    inline const double& yy() const {return mat[4];}
    inline double& yy() {return mat[4];}
    
    inline const double& yz() const {return mat[5];}
    inline double& yz() {return mat[5];}

    inline const double& zx() const {return mat[6];}
    inline double& zx() {return mat[6];}
    
    inline const double& zy() const {return mat[7];}
    inline double& zy() {return mat[7];}
    
    inline const double& zz() const {return mat[8];}
    inline double& zz() {return mat[8];}

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
    inline double tr() const {return mat[0] + mat[4] + mat[8];}

    // Get the invariants
    inline double I1() const {return tr();}

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
    // Contraction
    // return mat_ij source_ij
    // ------------------------------------------------------------------------
    double MatContraction( const Tensor2_3D &source ) const;

    // ------------------------------------------------------------------------
    // Contraction with the transpose of the input
    // return mat_ij source_ji
    // ------------------------------------------------------------------------
    double MatTContraction( const Tensor2_3D &source ) const;

    // print the matrix
    void print(std::ostream& os = std::cout, const std::string& delimiter = "\t") const;

    // print the matrix in a row
    void print_in_row() const;

    // print the Voigt components in the order of xx yy zz yz xz xy
    void print_Voigt() const;

    // ------------------------------------------------------------------------
    // Eigen decomposition of the matrix 
    // M = eta1 v1 v1T + eta2 v2 v2T + eta3 v3 v3T. 
    // This function will
    // return 1 if the three eigenvalues are the same
    // return 2 if there are two identical eigenvalues, the most distinct one is
    //          eta_1 and eta_1's associated eigenvector is v1
    // return 3 if all three are distinct
    // Reference: CMAME 197 2008 4007-4015 by W.M. Scherzinger and C.R. Dohrmann.
    // Note: this function is applicable to symmetric tensors only. However, we 
    // do not perform symmetry check within this function.
    // ------------------------------------------------------------------------
    int eigen_decomp( double &eta1, double &eta2, double &eta3,
       Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const;

  private:
    std::array<double,9> mat;

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
  inline Tensor2_3D transpose( const Tensor2_3D &input )
  {
    return Tensor2_3D( input(0), input(3), input(6),
      input(1), input(4), input(7),
      input(2), input(5), input(8) );
  }

  // Return an identity matrix
  inline Tensor2_3D gen_id()
  {
    return Tensor2_3D( 1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0 );
  }

  // Return a zero matrix
  inline Tensor2_3D gen_zero()
  {
    return Tensor2_3D( 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 );
  }

  // Return a random matrix
  inline Tensor2_3D gen_rand(double left = -1.0, double right = 1.0) 
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(left, right);
    
    return Tensor2_3D( dis(gen), dis(gen), dis(gen),
        dis(gen), dis(gen), dis(gen),
        dis(gen), dis(gen), dis(gen) );
    }

  // Return the exponential of the input matrix
  // exp(X) = sum_{k=0}^{infty} 1/(k!) X^k
  // The algorithm can be found on page 749 of Computational Methods for Plasticity
  // written by E.A. De Souza Neto, D. Peric and D.R.J. Owen.
  Tensor2_3D exp( const Tensor2_3D &input );
}

#endif
