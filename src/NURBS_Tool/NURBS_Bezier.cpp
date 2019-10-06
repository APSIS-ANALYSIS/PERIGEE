#include <iostream>
#include <assert.h>
#include <stdlib.h>

#include "NURBS_Bezier.hpp"

NURBS_T::BezierElem::BezierElem(const int sDegree,
                                const int tDegree,
                                const int uDegree,
                                bool initFullMat)
{
  mSDegree = sDegree;
  mTDegree = tDegree;
  mUDegree = uDegree;

  mNumCPts = (mSDegree+1)*(mTDegree+1)*(mUDegree+1);
  
  mDim = 1 + (tDegree > 0 ? 1 : 0) + (uDegree > 0 ? 1 : 0);
  
  mInitFullMat = initFullMat;
  mCoefMat = mCoefMatS = mCoefMatT = mCoefMatU = NULL;
  initMatrix();
}

NURBS_T::BezierElem::~BezierElem()
{
  int ii;

  if( NULL != mCoefMatS )
  {
    for( ii = 0; ii < mSDegree+1; ii++ )
      delete [] mCoefMatS[ii];
    delete [] mCoefMatS;
  }
  if( NULL != mCoefMatT )
  {
    for( ii = 0; ii < mTDegree+1; ii++ )
      delete [] mCoefMatT[ii];
    delete [] mCoefMatT;
  }
  
  if( NULL != mCoefMatU )
  {
    for( ii = 0; ii < mUDegree+1; ii++ )
      delete [] mCoefMatU[ii];
    delete [] mCoefMatU;
  }

  if( NULL != mCoefMat )
  {
    for( ii = 0; ii < mNumCPts; ii++ )
      delete [] mCoefMat[ii];
    delete [] mCoefMat;
  }
}

void NURBS_T::BezierElem::initMatrix()
{
  if( mInitFullMat )
  {
    mCoefMat = new double*[mNumCPts];
    
    int ii, jj;
    for( ii = 0; ii < mNumCPts; ii++ )
    {
      mCoefMat[ii] = new double[mNumCPts];
      
      //Initialize values (diagonal matrix).
      for( jj = 0; jj < mNumCPts; jj++ )
        mCoefMat[ii][jj] = 0.0;
      mCoefMat[ii][ii] = 1.0;
    }
  }
  else 
  {
    mCoefMatS = new double*[mSDegree+1];
    
    int ii, jj;
    for( ii = 0; ii < mSDegree+1; ii++ )
    {
      mCoefMatS[ii] = new double[mSDegree+1];
      
      //Initialize values (diagonal matrix).
      for( jj = 0; jj < mSDegree+1; jj++ )
        mCoefMatS[ii][jj] = 0.0;
      mCoefMatS[ii][ii] = 1.0;
    }
    
    if( mDim > 1 )
    {
      mCoefMatT = new double*[mTDegree+1];
      for( ii = 0; ii < mTDegree+1; ii++ )
      {
        mCoefMatT[ii] = new double[mTDegree+1];
        
        //Initialize values (diagonal matrix).
        for( jj = 0; jj < mTDegree+1; jj++ )
          mCoefMatT[ii][jj] = 0.0;
        mCoefMatT[ii][ii] = 1.0;
      }
      if( mDim > 2 )
      {
        mCoefMatU = new double*[mUDegree+1];
        for( ii = 0; ii < mUDegree+1; ii++ )
        {
          mCoefMatU[ii] = new double[mUDegree+1];
          
          //Initialize values (diagonal matrix).
          for( jj = 0; jj < mUDegree+1; jj++ )
            mCoefMatU[ii][jj] = 0.0;
          mCoefMatU[ii][ii] = 1.0;
        }
      }
    }
  }
}

double NURBS_T::BezierElem::coefficient(const int row,
                                        const int col) const
{
  if( row >= mNumCPts )
  {
    std::cerr << "Error: Row index out of bounds in NURBS_T::BezierElem::coefficient()" << std::endl;
    exit(1);
  }
  if( col >= mNumCPts )
  {
    std::cerr << "Error: Column index out of bounds in NURBS_T::BezierElem::coefficient()" << std::endl;
    exit(1);
  }
  
  if( NULL == mCoefMat )
  {
    if( 1 == mDim )
    {
      return mCoefMatS[row][col];
    }
    else if( 2 == mDim )
    {
      int s_Row = row % (mSDegree+1);
      int t_Row = (row - s_Row) / (mSDegree+1);
      assert( 0 == (row - s_Row) % (mSDegree+1) );
      
      int s_Col = col % (mSDegree+1);
      int t_Col = (col - s_Col) / (mSDegree+1);
      assert( 0 == (col - s_Col) % (mSDegree+1) );
      
      return mCoefMatS[s_Row][s_Col]*mCoefMatT[t_Row][t_Col];
    }
    else 
    {
      //Level 0 -> Full matrix
      //Level 1 -> First outer product blocks
      //Level 2 -> Univariat submatrix blocks.
      int level_1_Row = row % ( (mSDegree+1) * (mTDegree+1) );
      
      int s_Row = level_1_Row % (mSDegree+1);
      int t_Row = ( level_1_Row - s_Row ) / (mSDegree+1);
      int u_Row = ( row - s_Row - ( t_Row * (mSDegree+1) ) ) / ( (mSDegree+1) * (mTDegree+1) );
      
      int level_1_Col = col % ( (mSDegree+1) * (mTDegree+1) );
      int s_Col = level_1_Col % (mSDegree+1);
      int t_Col = ( level_1_Col - s_Col ) / (mSDegree+1);
      int u_Col = ( col - s_Col - ( t_Col * (mSDegree+1) ) ) / ( (mSDegree+1) * (mTDegree+1) );
      
      return mCoefMatS[s_Row][s_Col]*mCoefMatT[t_Row][t_Col]*mCoefMatU[u_Row][u_Col];
    }
  }
  else
    return mCoefMat[row][col];
}

void NURBS_T::BezierElem::set_S_MatrixCoeff(const int row,
                                            const int col,
                                            const double setValue)
{
  if( row >= mSDegree+1 )
  {
    std::cerr << "Error: Row index out of bounds in NURBS_T::BezierElem::set_S_MatrixCoeff()" << std::endl;
    exit(1);
  }
  if( col >= mSDegree+1 )
  {
    std::cerr << "Error: Column index out of bounds in NURBS_T::BezierElem::set_S_MatrixCoeff()" << std::endl;
    exit(1);
  }   
  
  mCoefMatS[row][col] = setValue;
  
}
void NURBS_T::BezierElem::set_T_MatrixCoeff(const int row,
                                            const int col,
                                            const double setValue)
{
  if( row >= mTDegree+1 )
  {
    std::cerr << "Error: Row index out of bounds in NURBS_T::BezierElem::set_T_MatrixCoeff()" << std::endl;
    exit(1);
  }
  if( col >= mTDegree+1 )
  {
    std::cerr << "Error: Column index out of bounds in NURBS_T::BezierElem::set_T_MatricCoeff()" << std::endl;
    exit(1);
  }  
  
  mCoefMatT[row][col] = setValue;
}
void NURBS_T::BezierElem::set_U_MatrixCoeff(const int row,
                                            const int col,
                                            const double setValue)
{
  if( row >= mUDegree+1 )
  {
    std::cerr << "Error: Row index out of bounds in NURBS_T::BezierElem::set_U_MatrixCoeff()" << std::endl;
    exit(1);
  }
  if( col >= mUDegree+1 )
  {
    std::cerr << "Error: Column index out of bounds in NURBS_T::BezierElem::set_U_MatrixCoeff()" << std::endl;
    exit(1);
  }  
  
  mCoefMatU[row][col] = setValue;
}

void NURBS_T::BezierElem::coefficient(const int row,
                                      const int col,
                                      const double newVal)
{
  if( NULL == mCoefMat && 1 != mDim )
  {
    std::cerr << "ERROR: Cannot add coefficient to multi-variate element matrix in "
    << "NURBS_T::BezierElem::coefficient()." << std::endl;
    exit(1);
  }
  
  if( row >= mNumCPts )
  {
    std::cerr << "Error: Row index out of bounds in NURBS_T::BezierElem::coefficient()" << std::endl;
    exit(1);
  }
  if( col >= mNumCPts )
  {
    std::cerr << "Error: Column index out of bounds in NURBS_T::BezierElem::coefficient()" << std::endl;
    exit(1);
  }  
  
  if( NULL == mCoefMat )
    mCoefMatS[row][col] = newVal;
  else
    mCoefMat[row][col] = newVal;
}

//EOF
