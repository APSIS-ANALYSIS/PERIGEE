#ifndef SYMMTENSOR4_3D_HPP
#define SYMMTENSOR4_3D_HPP
// ============================================================================
// SymmTensor4_3D.hpp
//
// This is a symmetric rank-4 tensor class. Here symmetric means both the major
// and minor symmetry.
// There are 6+5+...+1=21 entries in this object, which are stored in an array. 
// The ordering follows the Voigt notation.
// 
// ten[0]  -> I=0, J=0 : 0000
// ten[1]  -> I=0, J=1 : 0011 1100
// ten[2]  -> I=0, J=2 : 0022 2200
// ten[3]  -> I=0, J=3 : 0012 0021 1200 2100
// ten[4]  -> I=0, J=4 : 0002 0020 0200 2000
// ten[5]  -> I=0, J=5 : 0001 0010 0100 1000
// 
// ten[6]  -> I=1, J=1 : 1111
// ten[7]  -> I=1, J=2 : 1122 2211
// ten[8]  -> I=1, J=3 : 1112 1121 1211 2111
// ten[9]  -> I=1, J=4 : 1102 1120 0211 2011
// ten[10] -> I=1, J=5 : 1101 1110 0111 1011 
//
// ten[11] -> I=2, J=2 : 2222
// ten[12] -> I=2, J=3 : 2212 2221 1222 2122
// ten[13] -> I=2, J=4 : 2202 2220 0222 2022
// ten[14] -> I=2, J=5 : 2201 2210 0122 1022
// 
// ten[15] -> I=3, J=3 : 1212 1221 2112 2121
// ten[16] -> I=3, J=4 : 1202 1220 2102 2120 0212 2012 0221 2021
// ten[17] -> I=3, J=5 : 1201 1210 2101 2110 0112 1012 0121 1021
// 
// ten[18] -> I=4, J=4 : 0202 0220 2002 2020
// ten[19] -> I=4, J=5 : 0201 0210 2001 2010 0102 1002 0120 1020
//
// ten[20] -> I=5, J=5 : 0101 0110 1001 1010
//
// Author: Ju Liu
// Date: Aug. 3rd 2023
// ============================================================================
#include "Tensor4_3D.hpp"

class SymmTensor4_3D
{
  public:
    // ------------------------------------------------------------------------
    // Default constructor:
    // It construct the tensor with components ijkl being 
    // 0.5 * (delta_ik delta_jl + delta_il delta_jk) = dA_ij / dA_kl
    // with A = A^T.
    // ------------------------------------------------------------------------
    SymmTensor4_3D();

    ~SymmTensor4_3D();

    // Assignment operator
    SymmTensor4_3D& operator= (const SymmTensor4_3D &source);

    // Parenthesis operator: access through single index with 0 <= index < 21
    double& operator()(const int &index) {return ten[index];}

    const double& operator()(const int &index) const {return ten[index];}

    bool is_identical(const Tensor4_3D &source, const double &tol = 1.0e-12) const;
    
    bool is_identical(const SymmTensor4_3D &source, const double &tol = 1.0e-12) const;

    void print() const;

    void print_in_mat() const;

    // Addition operator : return left + right
    friend SymmTensor4_3D operator+( const SymmTensor4_3D &left, const SymmTensor4_3D &right);

    // Minus operator : return left - right
    friend SymmTensor4_3D operator-( const SymmTensor4_3D &left, const SymmTensor4_3D &right);

    // Add the source tensor to the object
    SymmTensor4_3D& operator+=( const SymmTensor4_3D &source );

    // Minus the source tensor to the object
    SymmTensor4_3D& operator-=( const SymmTensor4_3D &source );

    // Scalar multiplication
    SymmTensor4_3D& operator*=( const double &val );

    void gen_rand();

    void gen_zero();
    
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
    // add an outer product with scaling factor:
    //            ten_ijkl += val * mmat_ij  * mmat_kl
    // This is often used in the evaluation of the stiffness tensor.
    // ------------------------------------------------------------------------
    void add_OutProduct( const double &val, const SymmMatrix_3x3 &mmat );

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
    // Notes: vec1 = vec3 and vec2 = vec4 to maintain the major symmetry
    // ------------------------------------------------------------------------
    void add_SymmOutProduct( const double &val, const Vector_3 &vec1,
        const Vector_3 &vec2, const Vector_3 &vec3, const Vector_3 &vec4 );

    // ------------------------------------------------------------------------
    // add a symmetric product with a scaling factor -- val:
    // ten_ijkl += val * [ 0.5 * (mleft_ik mright_jl + mleft_il mright_jk) ]
    // This is often used in the evaluation of the stiffness tensor.
    // E.G., partial C^{-1}_AB / partial C_CD
    //     = -0.5 (C^{-1}_AC C^{-1}_BD + C^{-1}_AD C^{-1}_{BC})
    //     = SymmProduct(-0.5, invC, invC )
    // for invertible and symmetric 2nd-order tensor C.
    // Holzapfel book, p. 254
    // Notes: taking same mleft and mright maintains major symmetry, another
    // fact that mleft and mright are all symmetric maintains minor symmetry 
    // ------------------------------------------------------------------------
    void add_SymmProduct( const double &val, const SymmMatrix_3x3 &mleft,
        const SymmMatrix_3x3 &mright );

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
    // Notes: mleft and mright being all symmetric maintains minor symmetry
    // ------------------------------------------------------------------------
    void add_SymmOutProduct( const double &val, const SymmMatrix_3x3 &mleft,
        const SymmMatrix_3x3 &mright );

    // ------------------------------------------------------------------------
    // transform the natural indices of forth-order symmetric tensor to Voigt 
    // notation, minor symmetry requires ij / ji: 3x3 -> 6, kl / lk: 3x3 -> 6,
    // major symmetry requires ij_kl / ij_lk / ji_kl / ji_lk: 6x6 -> 21,
    // for more information, check the diagram above.
    // ------------------------------------------------------------------------
    int Voigt_notation( const int &ii, const int &jj, const int &kk, const int &ll ) const;

    // ------------------------------------------------------------------------
    // transfrom the natural indices of 3x3 symmetric matrix to Voigt notation,
    // so as to say ij / ji: 3x3 -> 6
    // ------------------------------------------------------------------------ 
    int Voigt_notation( const int &ii, const int &jj ) const;

  private:
    double ten[21];
};
#endif
