#ifndef MATRIX_SEQDENSE_PETSC32P7_HPP
#define MATRIX_SEQDENSE_PETSC32P7_HPP
// ==================================================================
// Matrix_SeqDense_Petsc.hpp
// Description:
// PETSc seq dense matrix as a derivde class of Matrix_SeqDense
// The Petsc lib should be in 3.2p7 version. 
// For other versions of PETSc, I need to modify some function call
// in this class.
//
// Date:
// Oct. 24 2013
// ==================================================================
#include "Matrix_double.hpp"
#include "Matrix_SeqDense.hpp"
#include "petscksp.h"

class Matrix_SeqDense_Petsc32p7 : public Matrix_SeqDense
{
  public:
    Matrix_SeqDense_Petsc32p7( const class Matrix_double * const &in_mat );
    virtual ~Matrix_SeqDense_Petsc32p7();

    virtual void print() const;

    virtual void solve( const vector<double> &b, vector<double> &x) const;

  private:
    Mat A;
    KSP ksp;
    PC pc;
};
#endif
