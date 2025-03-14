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
    constexpr SymmTensor2_3D() : mat{{1.0, 1.0, 1.0, 0.0, 0.0, 0.0}} {}

    // Copy constructor
    SymmTensor2_3D( const SymmTensor2_3D &source ) : mat(source.mat) {}

    // Constructor by six numbers in Voigt numbering
    constexpr SymmTensor2_3D( const double &m0, const double &m1, 
        const double &m2, const double &m3, const double &m4, const double &m5 ) 
      : mat{{m0, m1, m2, m3, m4, m5}} {}

    // Destructor
    ~SymmTensor2_3D() = default;

    // Assignment operator
    SymmTensor2_3D& operator= (const SymmTensor2_3D &source);

    inline Tensor2_3D full() const
    {
      return Tensor2_3D( mat[0], mat[5], mat[4], mat[5], mat[1], mat[3], mat[4], mat[3], mat[2] );
    }

    std::vector<double> to_std_vector() const 
    { return std::vector<double>(std::begin(mat), std::end(mat)); }

    std::array<double,6> to_std_array() const {return mat;}

    // Parenthesis operator. It allows accessing and assigning the matrix entries.
    inline double& operator()(const int &index) {return mat[index];}

    inline const double& operator()(const int &index) const {return mat[index];}

    // Get-functions that access components directly
    inline const double& xx() const {return mat[0];}
    inline double& xx() {return mat[0];}

    inline const double& xy() const {return mat[5];}
    inline double& xy() {return mat[5];}

    inline const double& xz() const {return mat[4];}
    inline double& xz() {return mat[4];}

    inline const double& yx() const {return mat[5];}
    inline double& yx() {return mat[5];}

    inline const double& yy() const {return mat[1];}
    inline double& yy() {return mat[1];}

    inline const double& yz() const {return mat[3];}
    inline double& yz() {return mat[3];}

    inline const double& zx() const {return mat[4];}
    inline double& zx() {return mat[4];}

    inline const double& zy() const {return mat[3];}
    inline double& zy() {return mat[3];}

    inline const double& zz() const {return mat[2];}
    inline double& zz() {return mat[2];}
    
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

    // unary minus operator
    SymmTensor2_3D operator-() const;

    // Return true if the input matrix is identical to the mat
    bool is_identical( const SymmTensor2_3D &source, const double &tol = 1.0e-12 ) const;

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
    inline double tr() const {return mat[0] + mat[1] + mat[2];}

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

    // Matrix rotation
    // Let Q be a rotation matrix, the matrix gets updated by
    // Q^T M Q = Q_ki M_kl Q_lj = output_matrix_ij
    void MatRot( const Tensor2_3D &Q );

    // Push-forward for contravariant tensor (like stress)
    // defined as F (object) Ft, see Holzapfel book p. 83
    void push_forward_stress( const Tensor2_3D &F );

    // Pull-back for contravariant tensor (like stress)
    // defined as invF (object) invFt, see Holzapfel book (2.87) 
    // in p. 83.
    void pull_back_stress( const Tensor2_3D &invF );    

    // Contraction
    // return mat_ij source_ij
    double MatContraction( const Tensor2_3D &source ) const;
    
    double MatContraction( const SymmTensor2_3D &source ) const;

    // Contraction
    // return source_ij source_ij
    double MatContraction() const;

    // Obtain the Voigt notation for a regular matrix index
    // 0 <= index < 9
    inline int Voigt_notation( const int &index ) const
    {return VoigtMap[index];}

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
    std::array<double,6> mat;

    // Define the Voigt map
    static constexpr std::array<int, 9> VoigtMap {{ 0, 5, 4, 5, 1, 3, 4, 3, 2 }};

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
  // Generate an identity tensor
  inline SymmTensor2_3D gen_id()
  { return SymmTensor2_3D(1.0, 1.0, 1.0, 0.0, 0.0, 0.0); }

  // Generate a zero tensor
  inline SymmTensor2_3D gen_zero()
  { return SymmTensor2_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }

  // Set components a random value
  inline SymmTensor2_3D gen_rand(const double &left = -1.0, const double &right = 1.0)
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(left, right);
    return SymmTensor2_3D( dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen) );
  }

  // Generate a dyad with a unit vector
  // Note: we do not check the unit length of input
  inline SymmTensor2_3D gen_dyad(const Vector_3 &input)
  {
    return SymmTensor2_3D( input(0) * input(0), input(1) * input(1),
        input(2) * input(2), input(1) * input(2), input(0) * input(2),
        input(0) * input(1) );
  }

  // Return the inverse of the input matrix
  SymmTensor2_3D inverse( const SymmTensor2_3D &input );

  // Generate the right Cauchy-Green tensor C = F^T F
  SymmTensor2_3D gen_right_Cauchy_Green( const Tensor2_3D &input );

  // Generate the left Cauchy-Green tensor b = F F^T
  SymmTensor2_3D gen_left_Cauchy_Green( const Tensor2_3D &input );

  // Convert a regular matrix to its symmetric part
  // output = 0.5 x ( source + source_transpose )
  SymmTensor2_3D gen_symm_part( const Tensor2_3D &input );

  // Apply the projector P := I - 1/3 invC x C on a symmetric tensor to obtain
  // its Dev part in the Lagrangian setting
  SymmTensor2_3D gen_DEV_part( const SymmTensor2_3D &input, const SymmTensor2_3D &CC );
}

#endif
