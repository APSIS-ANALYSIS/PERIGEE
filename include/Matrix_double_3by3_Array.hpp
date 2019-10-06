#ifndef MATRIX_DOUBLE_3BY3_ARRAY_HPP
#define MATRIX_DOUBLE_3BY3_ARRAY_HPP
// ============================================================================
// Matrix_double_3by3_Array.hpp
// This is a 3-by-3 matrix class that can calculate LU factorization of the
// dense matrix. The components are stored in a 1-D array.
//
// The array that stores the matrix is mat[9]. Logically, the matrix is 
//                    
//                     mat[0], mat[1], mat[2]
//                     mat[3], mat[4], mat[5]
//                     mat[6], mat[7], mat[8]
// 
// The p[3] array is a pointer for pivoting.
//
// Author: Ju Liu
// Date: Oct. 1st 2015
// ============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

class Matrix_double_3by3_Array
{
  public:
    // Defalt constructor: an identity matrix
    Matrix_double_3by3_Array();


    // Copy an input array and make it into a 3by3 matrix
    // Make sure in_mat has 9 components
    Matrix_double_3by3_Array(const double * const &in_mat);


    // Explicitly define the matrix components
    Matrix_double_3by3_Array( const double &a11, const double &a12, 
        const double &a13, const double &a21, const double &a22, 
        const double &a23, const double &a31, const double &a32, 
        const double &a33 );


    // Copy constructor
    Matrix_double_3by3_Array( const Matrix_double_3by3_Array &other );


    // Destructor
    ~Matrix_double_3by3_Array();


    // Assignment operator
    Matrix_double_3by3_Array& operator= (const Matrix_double_3by3_Array &input);


    // Parenthesis operator. It allows both access matrix entries as well as
    // assigning values to the entry.
    double& operator()(const int &index) {return mat[index];}
    
    const double& operator()(const int &index) const {return mat[index];}
    
    
    // Generate an identity matrix. Erase all previous values and reset p and
    // invm to default values.
    void gen_id();


    // Generate a matrix with random entries
    // All previous values are erased and p & invm are reset to default.
    void gen_rand();


    // Generate a Hilbert matrix
    // all previous values are earsed and p & invm are reset to default.
    void gen_hilb();


    // Perform LU factorization of the matrix object and store the L & U
    // matrix using the mat object.
    // Note: mat is changed after this call.
    void LU_fac();


    // Perform LU solve for the 3 mat x = b equations.
    // LU_fac() has to be called first.
    void LU_solve(const double * const &b, double * const &x ) const;


    // Perofrm LU solve for the 3 mat x = b equations
    // LU_fac() has to be called first.
    void LU_solve(const double &b1, const double &b2, const double &b3,
        double &x1, double &x2, double &x3) const;


    // Transpose operation for the matrix 
    void transpose();

    
    // Inverse of the matrix (based on cofactor). The p and invm are not
    // updated.
    void inverse();


    // determinant
    double det() const;


    // Vector multiplication y = Ax
    // make sure the x y vector has length 3.
    void VecMult( const double * const &x, double * const &y ) const; 

    
    // Matrix multiplication
    void MatMult( const Matrix_double_3by3_Array &mleft, 
        const Matrix_double_3by3_Array &mright );


    // print mat in matrix format
    void print() const;


    // print mat in matrix format and p & invm0, invm1, invm2
    void print_full() const;

  private:
    // Matrix entries
    double mat[9];
    
    // Pivoting flag. It is generated in LU-factorization and is used for
    // LU_solve.
    int p[3];
    
    // Inverse of the diagonal entries. It is generated in LU_fac and is 
    // used for LU_solve.
    double invm0, invm1, invm2;
};


#endif
