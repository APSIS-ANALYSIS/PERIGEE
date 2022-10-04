#ifndef TENSOR4_3D_HPP
#define TENSOR4_3D_HPP
// ============================================================================
// Tensor4_3D.hpp
// This is a 4th-order tensor in 3D. There are 3^4 = 81 double entries in this 
// object, which are stored in an array. The indices are arranged in the 
// following way:
//           ten[27i + 9j + 3k + l] = t_{ijkl}
// This mapping is inspired from the Matrix_3x3, where a matrix is mapped to an
// array as
//           mat[3i + j] = A_{ij}.
//
// The designed purpose is to ease the handling of stiffness tensor in solid 
// mechanics.
//
// Author: Ju Liu
// Date: July 3rd 2016
// ============================================================================
#include "Matrix_3x3.hpp"

class Tensor4_3D
{
  public:
    // ------------------------------------------------------------------------
    // Default constructor
    // Assuming delta_ij gives Konecker delta. The default 4th-order tensor
    // gives delta_ik delta_jl.
    // The intuition is that the A_{ijkl} b_{kl} = b_{ij}
    // This definition also comes from Holzapfel book, pp 23.
    // ------------------------------------------------------------------------
    Tensor4_3D();

    // Copy constructor
    Tensor4_3D( const Tensor4_3D &source );

    // Destructor
    ~Tensor4_3D();
    
    // Parenthesis operator: access through single index
    double& operator()(const int &index) {return ten[index];}

    const double& operator()(const int &index) const {return ten[index];}

    // Parenthesis operator: access through ii jj kk ll component index
    double& operator()(const int &ii, const int &jj, const int &kk, const int &ll)
    {return ten[27 * ii + 9 * jj + 3 * kk + ll];}

    const double& operator()(const int &ii, const int &jj, const int &kk, 
        const int &ll) const {return ten[27 * ii + 9 * jj + 3 * kk + ll];}

    bool is_identical(const Tensor4_3D &source, const double &tol = 1.0e-12) const;

    void print() const;

    // ------------------------------------------------------------------------
    // print the fourth-order tensor in the following matrix form:
    //
    //   CC_ijkl = [           l0        l1        l2
    //                  j0  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //             i0   j1  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //                  j2  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //
    //                  j0  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //             i1   j1  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //                  j2  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //
    //                  j0  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //             i2   j1  k0 k1 k2  k0 k1 k2  k0 k1 k2
    //                  j2  k0 k1 k2  k0 k1 k2  k0 k1 k2  ]
    // ------------------------------------------------------------------------
    void print_in_mat() const;

    // Copy operator
    void copy( const Tensor4_3D &source );

    // ------------------------------------------------------------------------
    // Generate basic 4th-order tensors.
    // ------------------------------------------------------------------------
    // Generate delta_ik delta_jl = dA_ij / dA_kl
    void gen_id();

    // ------------------------------------------------------------------------
    // Generate 0.5 * (delta_ik delta_jl + delta_il delta_jk) = dA_ij / dA_kl
    // with A = A^T.
    // Note: this is the derivative for symmetric 2nd-order tensor. In
    // principle, the derivative for symmetric tensor is nonunique, since the
    // derivative is acting on a symmetric tensor for the linearization and adding
    // a skew-symmetric tensor will not changing the effect. Hence, we define
    // the symmetric part of the 4th-order tensor be the derivative for the
    // 2nd-order tensor.
    // ------------------------------------------------------------------------
    void gen_symm_id();

    // ------------------------------------------------------------------------
    // Generate devaitoric projector
    // P_dev = Id4 - id2 corss id2 = delta_ik delta_jl - 1/3 delta_ij delta_kl
    // Holzapfel book, p. 24.
    // ------------------------------------------------------------------------
    void gen_proj_dev();

    // ------------------------------------------------------------------------
    // Generate Projector P = SymmId4 - 1/3 invC x C
    // P_IJKL = SymmID_IJKL - 1/3 invC_IJ C_KL
    // C is assumed to be the right Cauchy-Green tensor
    // invC is the inverse of C
    // see Holzapfel book p.229 eqn. (6.84).
    // ------------------------------------------------------------------------
    void gen_P( const Matrix_3x3 &C, const Matrix_3x3 &invC );

    // ------------------------------------------------------------------------
    // Generate Projector Ptilde = invC O invC - 1/3 invC x invC
    // invC is assumed to be the right Cauchy-Green tensor 
    // see Holzapfel book p. 255, eqn. (6.170).
    // ------------------------------------------------------------------------
    void gen_Ptilde( const Matrix_3x3 &invC );

    // ------------------------------------------------------------------------
    // generate a random 4th-order tensor (mainly used for debuggin)
    // ------------------------------------------------------------------------
    void gen_rand();

    // ------------------------------------------------------------------------
    // generate a zero 4th-order tensor
    // ------------------------------------------------------------------------
    void gen_zero();

    // ------------------------------------------------------------------------
    // Tensor algebraic manipulations: 
    // ------------------------------------------------------------------------
    // Scale the tensor by a scalar
    void scale( const double &val );

    // ten += input
    void PY( const Tensor4_3D &input );

    // ten += val * input
    void AXPY( const double &val, const Tensor4_3D &input );

    // ------------------------------------------------------------------------
    // add an outer product with scaling factor:
    //            ten_ijkl += val * mleft_ij  * mright_kl
    // This is often used in the evaluation of the stiffness tensor.
    // ------------------------------------------------------------------------
    void add_OutProduct( const double &val, const Matrix_3x3 &mleft,
        const Matrix_3x3 &mright );

    // ------------------------------------------------------------------------
    // add a tensor product of 4 vectors which is formed in the following way,
    //   vec1_i x vec2_j x vec3_k x vec4_l
    // This is typically used in the generation of the elasticity tensor for
    // stretch-based models, such as the Ogden model.
    // See, Holzapfel book p. 257. 
    // ------------------------------------------------------------------------
    void add_OutProduct( const double &val, const Vector_3 &vec1, const Vector_3 &vec2,
        const Vector_3 &vec3, const Vector_3 &vec4 );

    // ------------------------------------------------------------------------
    // add a symmetric product with a scaling factor -- val:
    // ten_ijkl += val * [ 0.5 * (mleft_ik mright_jl + mleft_il mright_jk) ]
    // This is often used in the evaluation of the stiffness tensor.
    // E.G., partial C^{-1}_AB / partial C_CD 
    //     = -0.5 (C^{-1}_AC C^{-1}_BD + C^{-1}_AD C^{-1}_{BC})
    //     = SymmProduct(-0.5, invC, invC )
    // for invertible and symmetric 2nd-order tensor C.
    // Holzapfel book, p. 254
    // ------------------------------------------------------------------------
    void add_SymmProduct( const double &val, const Matrix_3x3 &mleft,
        const Matrix_3x3 &mright );

    // ------------------------------------------------------------------------
    // Matrix update for the tensor:
    // ------------------------------------------------------------------------
    //    ten[IJKL] = A[iI] (or A[jJ], A[kK], A[lL]) ten[IJKL]
    // Einstein notation applied for the above relation.
    // This is mainly designed for the push-forward/pull-back operator for the
    // tensor. A is often the deformation gradient.
    // MatMult_1 : ten[iJKL] = A[iI] ten[IJKL]
    void MatMult_1( const Matrix_3x3 &source );

    // MatMult_2 : ten[IjKL] = A[jJ] ten[IJKL]
    void MatMult_2( const Matrix_3x3 &source ); 

    // MatMult_3 : ten[IJkL] = A[kK] ten[IJKL]
    void MatMult_3( const Matrix_3x3 &source );

    // MatMult_4 : ten[IJKl] = A[lL] ten[IJKL]
    void MatMult_4( const Matrix_3x3 &source );

    // ------------------------------------------------------------------------
    // Contraction with a 2nd-order tensor
    // ------------------------------------------------------------------------
    // Left contraction: A_ij ten_ijkl = B_kl
    void LeftContraction( const Matrix_3x3 &source, Matrix_3x3 &out ) const;

    // Right contraction: ten_ijkl A_kl = B_ij 
    void RightContraction( const Matrix_3x3 &source, Matrix_3x3 &out ) const;

    // Left & Right contraction A_ij ten_ijkl B_kl
    double LnRContraction( const Matrix_3x3 &Left, const Matrix_3x3 &Right ) const;

    // Tensor contraction in_ijkl ten_ijkl
    double Ten4Contraction( const Tensor4_3D &input ) const;

    // ------------------------------------------------------------------------
    // Tensor multiplication
    // ------------------------------------------------------------------------
    // ten = tleft * tright, in indicial form,
    // ten_IJKL = tleft_IJMN tright_MNKL
    void TenMult( const Tensor4_3D &tleft, const Tensor4_3D &tright );

    // Tensor Right Multiplication
    // ten_IJKL = ten_IJMN tright_MNKL
    void TenRMult( const Tensor4_3D &tright );

    // Tensor Left Multiplication
    // ten_IJKL = tleft_IJMN ten_MNKL
    void TenLMult( const Tensor4_3D &tleft );

    // Tensor Left and Right Multiplication modification
    // ten_IJKL = tleft_IJMN ten_MNST tright_STKL
    void TenLRMult( const Tensor4_3D &tleft, const Tensor4_3D &tright );

    // Tensor Left and Right Multiplication modification with the same
    // tensor P. See Holzapfel p. 255, the first term in eqn. (6.168) for
    // an example of this function.
    // ten_IJKL = P_IJMN ten_MNST P_KLST
    void TenPMult( const Tensor4_3D &P );

    // Check the minor / major symmetries of the tensor
    // major symmetry: CC_ijkl = CC_klij
    bool is_major_sym( const double &tol = 1.0e-12 ) const;
    
    // minor symmetry: CC_ijkl = CC_ijlk and CC_ijkl = CC_jikl
    bool is_minor_sym( const double &tol = 1.0e-12 ) const;

  private:
    double ten[81];
};

#endif
