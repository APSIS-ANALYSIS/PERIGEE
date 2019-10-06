#ifndef MATRIX_DOUBLE_3BY3_HPP
#define MATRIX_DOUBLE_3BY3_HPP
// ========================================================
// Matrix_double_3by3.hpp
// Description:
// Linear Algebra for 3by3 matrices.
//
// Date:
// Oct. 24 2013
// ========================================================
#include <cmath>
#include "Matrix_double.hpp"

class Matrix_double_3by3 : public Matrix_double
{
  public:
    Matrix_double_3by3();
    Matrix_double_3by3(const vector<double> &input);

    virtual ~Matrix_double_3by3();

    void Transpose();

    virtual double Det() const;

    virtual void Invert();
};
#endif
