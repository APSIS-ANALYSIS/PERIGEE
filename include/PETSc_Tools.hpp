#ifndef PETSC_TOOLS_HPP
#define PETSC_TOOLS_HPP
// ==================================================================
// PETSc_Tools.hpp
//
// PETSc_Tools is a namespace containing a suite of tools for
// PETSc objects.
//
// Author: Ju Liu
// Date: Jan 19 2018
// ==================================================================
#include "Vec_Tools.hpp"

namespace PETSc_T
{
  // ----------------------------------------------------------------
  // Display PETSc object
  // ----------------------------------------------------------------
  // Print the sparse matrix
  inline void print( const Mat &K )
  {MatView(K, PETSC_VIEWER_STDOUT_WORLD);}

  // Print the vector
  inline void print( const Vec &vv )
  {VecView(vv, PETSC_VIEWER_STDOUT_WORLD);}

  // Print the index set
  inline void print( const IS &is )
  {ISView(is, PETSC_VIEWER_STDOUT_WORLD);}

  // Draw nonzero structure graphically
  inline void Mat_DrawNZ( const Mat &K )
  {MatView(K, PETSC_VIEWER_DRAW_WORLD);}

  // Display the info object for the Mat K in the local portion of 
  // cpu rank given.
  void MatInfo_Display_local( const Mat &K, const PetscMPIInt &rank );

  // Display the info object for the Mat K for global sum.
  void MatInfo_Display_global( const Mat &K );

  // Calculate the dnz and onz number for a matrix.
  // dnz and onz have length of the local number of rows belonging to the
  // matrix.
  void Get_dnz_onz( const Mat &K, std::vector<int> &dnz, 
      std::vector<int> &onz );

  // ----------------------------------------------------------------
  // Set options for a matrix.
  // ----------------------------------------------------------------
  // Fix the nonzero structure of the matrix
  inline void Fix_nonzero_str( const Mat &K )
  {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}

  // New allocation forbidden. Adding or inserting in new locations
  // will generate an error message.
  inline void Fix_nonzero_err_str( const Mat &K )
  {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}

  // Ignore new allocation. Adding or inserting in a new allocation will
  // NOT generate error message.
  inline void Release_nonzero_err_str( const Mat &K )
  {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}

  // Keep the nonzero structure of the matrix
  inline void Keep_nonzero_pattern( const Mat &K )
  {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);}

  // ----------------------------------------------------------------
  // Operate on matrix
  // ----------------------------------------------------------------
  // generate a sparse identity matrix K with local length lrow.
  // User is responsible to free Mat by calling MatDestroy(&K).
  void MatCreateId(Mat &K, const PetscInt &lrow);

  // Output the diagonal value of K to the vector diag
  inline void GetDiagonal(const Mat &K, Vec &diag)
  {MatGetDiagonal(K, diag);}

  // Perform K = LEFT K RIGHT where LEFT and RIGHT are two diagonal
  // matrices with entries stored as two vectors left and right.
  inline void DiagonalScale( const Mat &K, const Vec &left, const Vec &right)
  {MatDiagonalScale(K, left, right);}

  // Get the number of rows in the local portion for the matrix K.
  inline int GetNumLocalRow( const Mat &K )
  {
    int m;
    MatGetLocalSize(K, &m, NULL);
    return m;
  }

  // Get the number of columns in the local portion for the matrix K.
  inline int GetNumLocalColumn( const Mat &K )
  {
    int n;
    MatGetLocalSize(K, NULL, &n);
    return n;
  }

  // ----------------------------------------------------------------
  // Operate on vector
  // ----------------------------------------------------------------
  // x[ii] = a[ii] * x[ii] 
  inline void DiagonalScale( const Vec &a, Vec &x)
  {VecPointwiseMult(x,a,x);}

  // Get the |diag|^{-0.5}. If the entry value is small, correct the diag
  // value to be 1.0.
  void MinusSqrtVec(Vec &diag, const double &tol = 1.0e-15);

  // Get |diag|^{-1}. If the entry value is smaller than the tol value, 
  // assign the resulting value to be 1.0.
  void InvAbsVec( Vec &diag, const double &tol = 1.0e-15 ); 

  // return the value of the vector at a location.
  // It can ONLY get the values in the local portion. Use with Caution!
  double GetValue( const Vec &a, const int ii );

  // Return sqrt(a dot a)
  inline double Get2Norm( const Vec &a )
  {
    double val;
    VecNorm(a, NORM_2, &val);
    return val;
  }

  // Write a vector to disk with file_name
  void WriteBinary( const Vec &a, const char * const &file_name );
}

#endif
