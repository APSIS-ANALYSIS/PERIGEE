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
#include <vector>
#include "Sys_Tools.hpp"
#include "petsc.h"

namespace PETSc_T
{
  std::string get_version();

  // ----------------------------------------------------------------
  // Display PETSc object
  // ----------------------------------------------------------------
  // Print the sparse matrix
  inline void print( const Mat &K )
  {MatView(K, PETSC_VIEWER_STDOUT_WORLD);}

  // Print the vector
  inline void print( const Vec &vv )
  {VecView(vv, PETSC_VIEWER_STDOUT_WORLD);}

  // Save the sparse matrix to a file
  inline void write_to_file( const Mat &K, const std::string &filename ) 
  { 
    PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD; 
    PetscViewerSetType(viewer, PETSCVIEWERASCII); 
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); 
    PetscViewerFileSetName(viewer, filename.c_str()); 
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); 
    MatView(K, viewer); 
  } 

  // Save the vector to a file
  inline void write_to_file( const Vec &vv, const std::string &filename ) 
  { 
    PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD; 
    PetscViewerSetType(viewer, PETSCVIEWERASCII); 
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); 
    PetscViewerFileSetName(viewer, filename.c_str()); 
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); 
    VecView(vv, viewer); 
  } 

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
  void Get_dnz_onz(const Mat &K, std::vector<int> &dnz, std::vector<int> &onz);

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

  // Get the vector's local copy's length 
  inline int GetLocalSize( const Vec &vv )
  {
    PetscInt vvsize; VecGetLocalSize( vv, &vvsize ); 
    return static_cast<int>(vvsize);
  }
  
  // Get the ghost vector's *local plus ghost* length.
  int GetLocalGhostSize( const Vec &vv );

  // Get the *ghost* vector's local copy in vv_array. One is responsible for
  // allocating the vv_array before calling this function and freeing the memory
  // after usage. The size of the vv_array can be obtained from
  // PETSc_T::GetLocalGhostSize.
  void GetLocalArray( const Vec &vv, double * const &vv_array );

  std::vector<double> GetLocalArray( const Vec &vv );

  // Scatter from a global vector
  void Scatter( const Vec &gg, PetscInt * &idx_from, const int &length, double * &values );
  // Warning:
  // VecScatterCreate can only be simultaneously run with all the CPUs that 
  // confirmed by PetscInitialize.
  // In other words, one should not use this function with only a part of CPUs,
  // this behavior will result in unknown failure with no error massege.
  // !!!!!!!!!!  Bad Code Example  !!!!!!!!!!
  //   (run with 5 CPUs)
  //   ...
  //   const int rank = SYS_T::get_MPI_rank();
  //   for(int ii=0; ii<rank+1; ++ii)
  //   {
  //     PETSc_T::Scatter(something...);
  //     some_functions_to_do_something();
  //   }
  //   ...
  // Since rank 0 have only one chance to run PETSc_T::Scatter, but the others
  // attempt to run it again, this code will fail.
  // ---------- Solution ----------
  //   ...
  //   const int rank = SYS_T::get_MPI_rank();
  //   for(int ii=0; ii<rank+1; ++ii)
  //   {
  //     PETSc_T::Scatter(something...);
  //     some_functions_to_do_something();
  //   }
  //   const int size = SYS_T::get_MPI_size();
  //   for(int ii=0; ii<size-rank-1; ++ii)
  //   {
  //     PETSc_T::Scatter(anything...);
  //     (Do nothing but just synchronize all the CPUs)
  //   }

  // Write a vector to disk with file_name
  void WriteBinary( const Vec &a, const std::string &file_name );

  // ==========================================================================
  // The followings are system function to monitor system memory usages
  // dynamically.
  // Note: Before calling the following four functions, the user need to call
  // PetscMemorySetGetMaximumUsage() immediately after PetscInitialize.
  // ==========================================================================
  // --------------------------------------------------------------------------
  // 1. Print the Maximum memory used for the program. The usage is
  //    reduced to CPU 0 by MPI_SUM.
  // --------------------------------------------------------------------------
  inline void print_MaxMemUsage()
  {
    PetscLogDouble memo = 0.0, memototal = 0.0;
    PetscMemoryGetMaximumUsage(&memo);
    MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    if( !SYS_T::get_MPI_rank() )
    {
      std::cout<<"\n Maximum Memeory usage : ";
      SYS_T::print_mem_size(memototal); std::cout<<"\n";
    }
  }

  // --------------------------------------------------------------------------
  // 2. Print the Current memory used for the program. The usage is
  //    reduced to CPU 0 by MPI_SUM.
  // --------------------------------------------------------------------------
  inline void print_CurMemUsage()
  {
    PetscLogDouble memo = 0.0, memototal = 0.0;
    PetscMemoryGetCurrentUsage(&memo);
    MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    if( !SYS_T::get_MPI_rank() )
    {
      std::cout<<"\n Current Memeory usage : ";
      SYS_T::print_mem_size(memototal); std::cout<<"\n";
    }
  }

  // --------------------------------------------------------------------------
  // 3. Print the Maximum space PETSc has allocated. This function
  //    should be used with the command line argument -malloc
  // --------------------------------------------------------------------------
  inline void print_MaxMallocUsage()
  {
    PetscLogDouble memo = 0.0, memototal = 0.0;
    PetscMallocGetMaximumUsage(&memo);
    MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    if( !SYS_T::get_MPI_rank())
    {
      std::cout<<"\n Maximum PETSc malloced : ";
      SYS_T::print_mem_size(memototal); std::cout<<"\n";
    }
  }

  // --------------------------------------------------------------------------
  // 4. Print the Current space PETSc has allocated. This function
  //    should be used with the command line argument -malloc
  // --------------------------------------------------------------------------
  inline void print_CurMallocUsage()
  {
    PetscLogDouble memo = 0.0, memototal = 0.0;
    PetscMallocGetCurrentUsage(&memo);
    MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    if( !SYS_T::get_MPI_rank() )
    {
      std::cout<<"\n Current PETSc malloced : ";
      SYS_T::print_mem_size(memototal); std::cout<<"\n";
    }
  }
}

#endif
