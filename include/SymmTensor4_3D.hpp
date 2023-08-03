#ifndef SYMMTENSOR4_3D_HPP
#define SYMMTENSOR4_3D_HPP
// ============================================================================
// SymmTensor4_3D.hpp
//
// This is a symmetric rank-4 tensor class. Here symmetric means both the major
// and minor symmetry.
// There are 6x6=36 entries in this object, which are stored in an array. 
// The ordering follows the Voigt notation.
//
// Author: Ju Liu
// Date: Aug. 3rd 2023
// ============================================================================

class SymmTensor4_3D
{
  public:
    // ------------------------------------------------------------------------
    // Default constructor:
    // It construct 
    // 0.5 * (delta_ik delta_jl + delta_il delta_jk) = dA_ij / dA_kl
    // with A = A^T.
    // ------------------------------------------------------------------------
    SymmTensor4_3D();

    ~SymmTensor4_3D();

    void gen_rand();

    void gen_zero();

  private:
    double ten[36];
};

#endif
