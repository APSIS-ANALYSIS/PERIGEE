#ifndef MATRIX_SEQDENSE_HPP
#define MATRIX_SEQDENSE_HPP
// ==================================================================
// Matrix_SeqDense.hpp
// Description:
// A template for manipulating sequential dense small matrix.
//
// Date:
// Oct. 24. 2013
// ==================================================================
#include <vector>

using namespace std;

class Matrix_SeqDense
{
  public:
    Matrix_SeqDense() {};
    virtual ~Matrix_SeqDense() {};

    virtual void print() const = 0;

    // solve the matrix eqaution: Ax = b
    virtual void solve(const vector<double> &b,
      vector<double> &x) const =0;
};
#endif
