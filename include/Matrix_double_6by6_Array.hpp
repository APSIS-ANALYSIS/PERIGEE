#ifndef MATRIX_DOUBLE_6BY6_ARRAY_HPP
#define MATRIX_DOUBLE_6BY6_ARRAY_HPP
// ==================================================================
// Matrix_double_6by6_Array.hpp
// This is a 6by6 matrix class that we used in calculating 2nd order
// derivatives of element basis functions.
//
// Date: Nov. 21 2013
// ==================================================================
#include <iostream>
#include <cmath>
#include <array>

class Matrix_double_6by6_Array
{
  public:
    Matrix_double_6by6_Array(const double * const &in_array);
    
    // This constructor is designed for the special structure when
    // solving for 2nd order derivatives.
    Matrix_double_6by6_Array(const double &a, const double &b,
        const double &c, const double &d, const double &e,
        const double &f, const double &g, const double &h,
        const double &i);
    
    ~Matrix_double_6by6_Array();

    void LU_fac();
    
    void LU_solve( const double * const &rhs, double * const &sol ) const;
    
    std::array<double, 6> LU_solve( const std::array<double, 6> &rhs ) const;

    void print() const;

  private:
    double Mat[6][6];
    int p[6];
    double invm0, invm1, invm2, invm3, invm4, invm5;
};

#endif
