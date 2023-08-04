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
    double ten[21];
};

#endif
