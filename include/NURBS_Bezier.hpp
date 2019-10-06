#ifndef NURBS_BEZIER_HPP
#define NURBS_BEZIER_HPP

#include <vector>

//#define STORE_UNIVARIATE_MATRICES

namespace NURBS_T
{
//! A Bezier Segment.    
    class BezierElem
    {
  public:
      BezierElem(const int sDegree,
                 const int tDegree = 0,
                 const int uDegree = 0,
                 bool initFullMat = true);
      ~BezierElem();

//! Stores map to global control points that influence this element.     
      std::vector<int> gPtsMap;
      
//! Stores map to local control points that define this element.      
      std::vector<int> lPtsMap;

//! Returns the coefficent value at (row,col) index.        
      double coefficient(const int row,
                         const int col) const;


     void set_S_MatrixCoeff(const int row,
                            const int col,
                            const double setValue);
     
      void set_T_MatrixCoeff(const int row,
                             const int col,
                             const double setValue);
      
      void set_U_MatrixCoeff(const int row,
                             const int col,
                             const double setValue);
      
//! Sets the coefficient value at (row,col) index.
      void coefficient(const int row,
                        const int col,
                        const double setValue);
 
//! Returns the polynomial degree of this Segment in the s direction.
      int sDegree()
          {return mSDegree;}
//! Returns the polynomial degree of this Segment in the t direction.
      int tDegree()
          {return mTDegree;}
//! Returns the polynomial degree of this Segment in the u direction.      
      int uDegree()
          {return mUDegree;}
      
  private:
      void initMatrix();
      bool mInitFullMat;
      
      int mDim;
      int mSDegree;
      int mTDegree;
      int mUDegree;

      int mNumCPts;
      double** mCoefMatS;
      double** mCoefMatT;
      double** mCoefMatU;
      double** mCoefMat;
    };
}
#endif
