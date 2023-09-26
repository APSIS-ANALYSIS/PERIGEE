#ifndef MATRIX_DOUBLE_6BY6_ARRAY_HPP
#define MATRIX_DOUBLE_6BY6_ARRAY_HPP
// ==================================================================
// Matrix_double_6by6_Array.hpp
// This is a 6by6 matrix class that we used in calculating 2nd order
// derivatives of element basis functions. The components are stored
// in a 1-D array.
//
// The array that stores the matrix is mat[36]. Logically, the matrix is 
//                    
//    mat[0],  mat[1],  mat[2],  mat[3],  mat[4],  mat[5]
//    mat[6],  mat[7],  mat[8],  mat[9],  mat[10], mat[11]
//    mat[12], mat[13], mat[14], mat[15], mat[16], mat[17]
//    mat[18], mat[19], mat[20], mat[21], mat[22], mat[23]
//    mat[24], mat[25], mat[26], mat[27], mat[28], mat[29]
//    mat[30], mat[31], mat[32], mat[33], mat[34], mat[35]
//
// Author: Ju Liu
// Date: Nov. 21 2013
// ==================================================================
#include <iostream>
#include <cmath>
#include <array>

class Matrix_double_6by6_Array
{
  public:
    // This constructor is designed for the special structure when
    // solving for 2nd order derivatives.
    Matrix_double_6by6_Array(const double &aa, const double &bb,
        const double &cc, const double &dd, const double &ee,
        const double &ff, const double &gg, const double &hh,
        const double &ii);
    
    ~Matrix_double_6by6_Array();

    void LU_fac();
    
    std::array<double, 6> LU_solve( const std::array<double, 6> &rhs ) const;

    void print() const;

  private:
    double Mat[36];
    int pp[6];
};

#endif
