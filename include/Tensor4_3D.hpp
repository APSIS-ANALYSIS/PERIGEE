#ifndef TENSOR4_3D_HPP
#define TENSOR4_3D_HPP
// ============================================================================
// Tensor4_3D.hpp
//
// This is a 4th-order tensor in 3D. There are 3^4 = 81 double entries in this 
// object, which are stored in an array. The indices are arranged in the 
// following way:
//           ten[27i + 9j + 3k + l] = t_{ijkl}
// This mapping is inspired from the Tensor2_3D, where a matrix is mapped to an
// array as
//           mat[3i + j] = A_{ij}.
//
// The design purpose is primarily to handle the stiffness tensor arising in 
// solid mechanics.
//
// Author: Ju Liu
// Date: July 3rd 2016
// ============================================================================
#include "SymmTensor2_3D.hpp"
#include "Math_Tools.hpp"

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
    Tensor4_3D() { gen_id(); }

    // Copy constructor
    Tensor4_3D( const Tensor4_3D &source ) : ten(source.ten) {}

    Tensor4_3D( const std::array<double, 81> &source ) : ten(source) {}

    // Destructor
    ~Tensor4_3D() = default;
   
    // Assignment operator
    Tensor4_3D& operator= (const Tensor4_3D &source);

    // Parenthesis operator: access through single index with 0 <= index < 81
    inline double& operator()(const int &index) {return ten[index];}

    inline const double& operator()(const int &index) const {return ten[index];}

    // Parenthesis operator: access through ii jj kk ll component index
    inline double& operator()(const int &ii, const int &jj, const int &kk, const int &ll)
    {return ten[27 * ii + 9 * jj + 3 * kk + ll];}

    inline const double& operator()(const int &ii, const int &jj, const int &kk, 
        const int &ll) const {return ten[27 * ii + 9 * jj + 3 * kk + ll];}

    inline std::vector<double> to_std_vector() const
    {return std::vector<double>(std::begin(ten), std::end(ten));}

    inline std::array<double,81> to_std_array() const {return ten;}

    bool is_identical(const Tensor4_3D &source, const double &tol = 1.0e-12) const;

    void print(std::ostream &os = std::cout, const std::string &delimiter = "\t") const;

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

    // Addition operator : return left + right
    friend Tensor4_3D operator+( const Tensor4_3D &left, const Tensor4_3D &right);

    // Minus operator : return left - right
    friend Tensor4_3D operator-( const Tensor4_3D &left, const Tensor4_3D &right);

    // Add the source tensor to the object
    Tensor4_3D& operator+=( const Tensor4_3D &source );

    // Minus the source tensor to the object
    Tensor4_3D& operator-=( const Tensor4_3D &source );

    // Scalar multiplication
    Tensor4_3D& operator*=( const double &val );

    // ------------------------------------------------------------------------
    // Generate basic 4th-order tensors.
    // ------------------------------------------------------------------------
    // Generate delta_ik delta_jl = dA_ij / dA_kl
    void gen_id();

    // ------------------------------------------------------------------------
    // Generate devaitoric projector
    // P_dev = Id4 - id2 corss id2 = delta_ik delta_jl - 1/3 delta_ij delta_kl
    // Holzapfel book, p. 24.
    // ------------------------------------------------------------------------
    void gen_proj_dev();

    // ------------------------------------------------------------------------
    // Tensor algebraic manipulations: 
    // ------------------------------------------------------------------------
    // Scale the tensor by a scalar
    void scale( const double &val );

    // ten += val * input
    void AXPY( const double &val, const Tensor4_3D &input );

    // Transpose a rank-four tensor ten_ijkl to be ten_klij, refered to the 
    // equation (1.160 ) in Holzapfel book, p.23.
    void transpose();  

    // ------------------------------------------------------------------------
    // add an outer product with scaling factor:
    //            ten_ijkl += val * mleft_ij  * mright_kl
    // This is often used in the evaluation of the stiffness tensor.
    // ------------------------------------------------------------------------
    void add_OutProduct( const double &val, const Tensor2_3D &mleft,
        const Tensor2_3D &mright );

    void add_OutProduct( const double &val, const SymmTensor2_3D &mleft,
        const SymmTensor2_3D &mright );

    // ------------------------------------------------------------------------
    // add a tensor product of 4 vectors which is formed in the following way,
    //         val x vec1[i] x vec2[j] x vec3[k] x vec4[l]
    // This is typically used in the generation of the elasticity tensor for
    // stretch-based models, such as the Ogden model.
    // See, Holzapfel book p. 257. 
    // ------------------------------------------------------------------------
    void add_OutProduct( const double &val,
        const Vector_3 &vec1, const Vector_3 &vec2,
        const Vector_3 &vec3, const Vector_3 &vec4 );

    // ------------------------------------------------------------------------
    // add a symmetric tensor product of 4 vectors which is defined the
    // following way,
    //     val x ( vex1[i] x vec2[j] x vec3[k] x vec4[l]
    //     + vex1[i] x vec2[j] x vec3[l] x vec4[k]
    //     + vex1[j] x vec2[i] x vec3[k] x vec4[l]
    //     + vex1[j] x vec2[i] x vec3[l] x vec4[k] ).
    // This function is typically called in the generation of the elasticity
    // tensor in the stretch-based models.
    // See, Holzapfel book p. 263, equation (6.196) for an example.
    // ------------------------------------------------------------------------
    void add_SymmOutProduct( const double &val, const Vector_3 &vec1, 
        const Vector_3 &vec2, const Vector_3 &vec3, const Vector_3 &vec4 );

    // ------------------------------------------------------------------------
    // add a symmetric product with a scaling factor -- val:
    // ten_ijkl += val * [ 0.5 * (mleft_ik mright_jl + mleft_il mright_jk) ]
    // This is often used in the evaluation of the stiffness tensor.
    // E.G., partial C^{-1}_AB / partial C_CD 
    //     = -0.5 (C^{-1}_AC C^{-1}_BD + C^{-1}_AD C^{-1}_{BC})
    //     = SymmProduct(-1.0, invC, invC )
    // for invertible and symmetric 2nd-order tensor C.
    // Holzapfel book, p. 254
    // ------------------------------------------------------------------------
    void add_SymmProduct( const double &val, const Tensor2_3D &mleft,
        const Tensor2_3D &mright );

    // ------------------------------------------------------------------------
    // add the out-product of two 2nd-order tensor in a symmetric fashion,
    // ten_ijkl += val * [  (mleft_ij mright_kl + mleft_kl mright_ij) ]
    // this is equivalent to and faster than
    //         add_OutProduct(val, mleft, mright); 
    //         add_OutProduct(val, mright, mleft);
    // This function can be used in the generation of the elasticity tensor. For
    // example, in Holzapfel book p. 261, the terms associalted with delta_2,
    // delta_3, and delta_5 in (6.193) can be generated with this function by
    //         add_SymmOutProduct(delta_2, I, C   );
    //         add_SymmOutProduct(delta_3, I, Cinv);
    //         add_SymmOutProduct(delta_5, C, Cinv);
    // Or, in the definition of C_iso of (6.168), this function can be used to
    // genereate the last term by
    //         add_SymmOutProduct(-2/3, Cinv, Siso);
    // ------------------------------------------------------------------------
    void add_SymmOutProduct( const double &val, const Tensor2_3D &mleft,
        const Tensor2_3D &mright );

    // ------------------------------------------------------------------------
    // Matrix update for the tensor:
    // ------------------------------------------------------------------------
    //    ten[IJKL] = A[iI] (or A[jJ], A[kK], A[lL]) ten[IJKL]
    // Einstein notation applied for the above relation.
    // This is mainly designed for the push-forward/pull-back operator for the
    // tensor. A is often the deformation gradient.
    // MatMult_1 : ten[iJKL] = A[iI] ten[IJKL]
    void MatMult_1( const Tensor2_3D &source );

    // MatMult_2 : ten[IjKL] = A[jJ] ten[IJKL]
    void MatMult_2( const Tensor2_3D &source ); 

    // MatMult_3 : ten[IJkL] = A[kK] ten[IJKL]
    void MatMult_3( const Tensor2_3D &source );

    // MatMult_4 : ten[IJKl] = A[lL] ten[IJKL]
    void MatMult_4( const Tensor2_3D &source );

    // ------------------------------------------------------------------------
    // Contraction with a 2nd-order tensor
    // ------------------------------------------------------------------------
    // Left contraction: A_ij ten_ijkl = B_kl
    Tensor2_3D LeftContraction( const Tensor2_3D &source ) const;

    // Right contraction: ten_ijkl A_kl = B_ij 
    Tensor2_3D RightContraction( const Tensor2_3D &source ) const;

    // Left & Right contraction A_ij ten_ijkl B_kl
    double LnRContraction( const Tensor2_3D &Left, const Tensor2_3D &Right ) const;

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

    // ------------------------------------------------------------------------
    // This function calculates AA_IJKL X_KL = B_IJ. In this function, the 
    // rank-four tensor AA needs to satisfy the minor symmetry, and the 
    // rank-two tensor B need to be symmetric. Then, using the Voigt notation,
    // we calculate a matrix problem AA_6x6 X_6 = B_6, and we then map the 
    // solution X_6 into the tensor X_KL. 
    // ------------------------------------------------------------------------
    Tensor2_3D solve( const Tensor2_3D &B ) const;
    
    SymmTensor2_3D solve( const SymmTensor2_3D &B ) const;

    // ------------------------------------------------------------------------
    // This function calculates AA_IJKL XX_KLMN = BB_IJMN. In this function,
    // the rank-four tensor BB will be divided into nine rank-two tensor B,
    // and the above function will be called to solve AA_IJKL X_KL = B_IJ.
    // The nine solutions X will be put together to form the solution XX.
    // ------------------------------------------------------------------------
    Tensor4_3D solve( const Tensor4_3D &BB ) const;

  private:
    std::array<double,81> ten;
};

// Right contraction: return ten_ijkl source_kl
Tensor2_3D operator*( const Tensor4_3D &ten, const Tensor2_3D &source );

// Left contraction: return source_ij ten_ijkl
Tensor2_3D operator*( const Tensor2_3D &source, const Tensor4_3D &ten );

// Tensor multiplication: return tleft_ijmn tright_mnkl
Tensor4_3D operator*( const Tensor4_3D &tleft, const Tensor4_3D &tright );

// Return scalar multiplication on the input tensor
Tensor4_3D operator*( const double &val, const Tensor4_3D &input );

namespace Ten4
{
  inline Tensor4_3D gen_zero()
  {
    std::array<double,81> out {};
    out.fill(0.0);
    return Tensor4_3D(out);
  }

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
  Tensor4_3D gen_symm_id();
    
  // ------------------------------------------------------------------------
  // generate a random 4th-order tensor (mainly used for debuggin)
  // ------------------------------------------------------------------------
  Tensor4_3D gen_rand(const double &left = -1.0, const double &right = 1.0);

  // ------------------------------------------------------------------------
  // Generate Projector P = SymmId4 - 1/3 invC x C
  // P_IJKL = SymmID_IJKL - 1/3 invC_IJ C_KL
  // C is assumed to be the right Cauchy-Green tensor
  // invC is the inverse of C
  // see Holzapfel book p.229 eqn. (6.84).
  // ------------------------------------------------------------------------
  Tensor4_3D gen_P( const Tensor2_3D &C, const Tensor2_3D &invC );

  Tensor4_3D gen_P( const Tensor2_3D &C );

  Tensor4_3D gen_P( const SymmTensor2_3D &C, const SymmTensor2_3D &invC );

  Tensor4_3D gen_P( const SymmTensor2_3D &C );

  // ------------------------------------------------------------------------
  // Generate Projector Pt = transpose of P = SymmId4 - 1/3 C x invC
  // P_IJKL = SymmID_IJKL - 1/3 C_IJ invC_KL
  // C is assumed to be the right Cauchy-Green tensor
  // invC is the inverse of C
  // see Holzapfel book p.229 eqn. (6.84).
  // ------------------------------------------------------------------------
  Tensor4_3D gen_Pt( const Tensor2_3D &C );
  
  Tensor4_3D gen_Pt( const SymmTensor2_3D &C );

  // ------------------------------------------------------------------------
  // Generate Projector Ptilde = invC O invC - 1/3 invC x invC
  // invC is assumed to be the right Cauchy-Green tensor 
  // see Holzapfel book p. 255, eqn. (6.170).
  // ------------------------------------------------------------------------
  Tensor4_3D gen_Ptilde( const Tensor2_3D &invC );

  // ------------------------------------------------------------------------
  // Generate a transpose of a rank-four tensor, refered to the equation (1.160) 
  // in Holzapfel book, p. 23.
  // ------------------------------------------------------------------------
  Tensor4_3D transpose( const Tensor4_3D &input );  
}

#endif
