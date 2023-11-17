#ifndef TENSOR4_2D_HPP
#define TENSOR4_2D_HPP
// ==================================================================
// Tensor4_2D.hpp
// This is a 4th-order tensor in 2D. There are 2^4 = 16 double entries
// in this object, which are stored in an array. The indices are 
// arranged in the following way:
//             ten[ 8i + 4j + 2k + l ] = t_{ijkl}
//
// The design is intended to ease the handling of stiffness tensor in
// 2D solid mechanics models.
//
// Author: Ju Liu
// Date: Sept. 12 2016
// ==================================================================
#include "Tensor2_2D.hpp" 

class Tensor4_2D
{
  public:
    // Default constructor
    // Assuming delta_ij is the Konecker delta. The default 4th-order
    // tensor is      t_ijkl =  delta_ik delta_jl.
    Tensor4_2D();

    // Copy constructor
    Tensor4_2D( const Tensor4_2D &source );

    // Destructor
    ~Tensor4_2D() = default;

    // Parenthesis operator: access through a single index
    double& operator()(const int &index) {return ten[index];}

    const double& operator()(const int &index) const {return ten[index];}

    // Parenthesis operator: access through ijkl component index
    double& operator()(const int &ii, const int &jj, const int &kk, const int &ll)
    {return ten[8 * ii + 4 * jj + 2 * kk + ll];}

    const double& operator()(const int &ii, const int &jj, const int &kk,
        const int &ll) const {return ten[8 * ii + 4 * jj + 2 * kk + ll];}

    void print() const;

    // Copy operator
    void copy( const Tensor4_2D &source );

    // Generate 4th-order tensor
    // t_ijkl = delta_ik delta_jl ( = dA_ij / dA_kl)
    void gen_id();

    // Generate 0.5 * (delta_ik delta_jl + delta_il delta_jk) = dA_ij / dA_kl
    // with A = A^T.
    // Note: this is the derivative for symmetric 2nd-order tensor. In
    // principle, the derivative for symmetric tensor is nonunique, since the
    // derivative is acting on a symmetric tensor for the linearization and adding
    // a skew-symmetric tensor will not changing the effect. Hence, we define
    // the symmetric part of the 4th-order tensor be the derivative for the
    // 2nd-order tensor.
    void gen_symm_id();

    // generate a random 4th-order tensor (mainly used for debugging)
    void gen_rand(const double &min = -1.0, const double &max = 1.0);

    // generate a zero 4th-order tensor
    void gen_zero();

    // Scale the tensor by a scalar
    void scale( const double &val );

    // ten += input
    void PY( const Tensor4_2D &input );

    // ten += val * input
    void AXPY( const double &val, const Tensor4_2D &input );

    // add an outer product with scaling factor:
    //            ten_ijkl += val * mleft_ij  * mright_kl
    // This is often used in the evaluation of the stiffness tensor.
    void add_OutProduct( const double &val, const Tensor2_2D &mleft,
        const Tensor2_2D &mright );

    // add a symmetric product with scaling factor:
    //  ten_ijkl += val * [ 0.5 * (mleft_ik mright_jl + mleft_il mright_jk) ]
    // This is often used in the evaluation of the stiffness tensor
    void add_SymmProduct( const double &val, const Tensor2_2D &mleft,
        const Tensor2_2D &mright );

    // Matrix update for the tensor:
    //    ten[IJKL] = A[iI] (or A[jJ], A[kK], A[lL]) ten[IJKL]
    // Einstein notation applied for the above relation.
    // This is mainly designed for the push-forward/pull-back operator for the
    // tensor. A is often the deformation gradient.
    // MatMult_1 : ten[iJKL] = A[iI] ten[IJKL]
    void MatMult_1( const Tensor2_2D &source );

    // MatMult_2 : ten[IjKL] = A[jJ] ten[IJKL]
    void MatMult_2( const Tensor2_2D &source ); 

    // MatMult_3 : ten[IJkL] = A[kK] ten[IJKL]
    void MatMult_3( const Tensor2_2D &source );

    // MatMult_4 : ten[IJKl] = A[lL] ten[IJKL]
    void MatMult_4( const Tensor2_2D &source );

    // Contraction with a 2nd-order tensor
    // Left contraction: A_ij ten_ijkl = B_kl
    void LeftContraction( const Tensor2_2D &source, Tensor2_2D &out ) const;

    // Right contraction: ten_ijkl A_kl = B_ij 
    void RightContraction( const Tensor2_2D &source, Tensor2_2D &out ) const;

    // Left & Right contraction A_ij ten_ijkl B_kl
    double LnRContraction( const Tensor2_2D &Left, 
        const Tensor2_2D &Right ) const;

    // Tensor contraction in_ijkl ten_ijkl
    double Ten4Contraction( const Tensor4_2D &input ) const;

  private:
    double ten[16];
};

#endif
