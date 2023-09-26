#ifndef MATRIX_DOUBLE_6BY6_ARRAY_HPP
#define MATRIX_DOUBLE_6BY6_ARRAY_HPP
// ==================================================================
// Matrix_double_6by6_Array.hpp
// This is a 6by6 matrix class that we used in calculating 2nd order
// derivatives of element basis functions.
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
    double Mat[6][6];
    int pp[6];
};

#endif
