#ifndef NURBS_TOOLS_HPP
#define NURBS_TOOLS_HPP
#include "NURBS_Bezier.hpp"
//!NUBS Tools namespace
namespace NURBS_T
{
  //   struct BezierElem;

  //! Determine knot span index
  /*!
    \param knotVec The knot vector
    \param cDegree The polynomial degree of the curve
    \param value The value for which the knot span will be determined
    \return The span index. Returns -1 if value is outside the knot vector end points.
    */
  int findSpan(const std::vector<double>& knotVec,
      const int cDegree,
      const double value);

  //! Project to homogeneous coordinates
  /*!
    Project the control points up one dimension to {x1*w1,y1*w1,z1*w1,w1,...}
    \param ctrlPts The points to project.
    \param spatialDim the spatial dimension of the non-homogeneous points (2 or 3)
    \return 0 if no errors, 1 otherwise.
    */
  int projectUp(std::vector<double>& ctrlPts,
      const int spatialDim);

  //! Project from homogeneous coordinates
  /*!
    Project the control points down one dimension to {x1,y1,z1,w1,...}
    \param ctrlPts The points to project.
    \param spatialDim the spatial dimension of the non-homogeneous points (2 or 3)
    \return 0 if no errors, 1 otherwise.
    */
  int projectDown(std::vector<double>& ctrlPts,
      const int spatialDim);

  //! Compute the nonvanishing basis functions
  /*!
    See Piegl and Tiller pg. 70
    \param knotVec The knot vector
    \param cDegree The degree of the curve
    \param u Parametric value at which to evaluate the basis functions
    must be within the end points of the knot vector
    \param N Return the cDegree+1 values of the basis functions.
    \return The result of findSpan() for u.
    */
  int basisFunctions(const std::vector<double>& knotVec,
      const int cDegree,
      const double u,
      std::vector<double>& N);

  //! Compute nonzero basis functions and their derivatives
  /*!
    See Piegl and Tiller Pg. 72
    \param knotVec The knot vector
    \param cDegree The degree of the curve
    \param u Parametric value at which to evaluate the basis functions.
    Must be within the end points of the knot vector.
    \param numDers The number of derivatives to calculate.  Note: if
    numDirs > cDegree will only calcualte cDegree derivatives.
    all others will be zero.
    \param ders Returns the computes values.  Values are orders as
    N0, dN0/dx, d^2N0/dx^2, ... , d^numDir N0/dx^numDirs, ...
    \return The result of findSpan() for u.
    */
  int basisFuncsAndDers(const std::vector<double>& knotVec,
      const int cDegree,
      const double u,
      const int numDers,
      std::vector<double>& ders);

  //! Compute the Greville abscissae for the knot vector.
  /*!
    \param knotVec The (open) knot vector
    \param cDegree The degree of the curve
    \param pts Returns the values of the Greville abscissae. 
    \return 0 if no errors, 1 otherwise
    */
  int computeGrevilleAbscissae(const std::vector<double>& knotVec,
      const int cDegree,
      std::vector<double>& pts);

  // Computes the point on a curve at fixed s values.
  /*!
    See Piegl and Tiller Pg. 82.  This function handles homogeneous
    coordinates so it works for Rational curves as well as B-Splines.
    \param knotVec The knot vector
    \param cDegree The degree of the curve
    \param Dim The dimension of the control points.
    \param ctrlPts The homogeneous control points for the surface
    \param s
    \param outPts A vector of the coordinates of the resulting point
    \return 0 if no errors, 1 otherwise
    */
  int curvePoint(const std::vector<double>& knotVec,
      const int cDegree,
      const int Dim,
      const std::vector<double>& ctrlPts,
      const double s,
      std::vector<double>& outPts);

  // Computes the point on a surface at fixed (s,t) values.
  /*!
    See Piegl and Tiller Pg. 103.  This function handles homogeneous
    coordinates so it works for Rational surfaces as well as B-Spline
    surfaces.
    \param sKnots The knots in the S direction
    \param sDegree The degree in the S direction
    \param tKnots The knots in the T direction
    \param tDegree The degree in the T direction
    \param Dim The dimension of the control points.
    \param ctrlPts The homogeneous control points for the surface
    \param s
    \param t
    \param outPts A vector of the coordinates of the resulting point
    \return 0 if no errors, 1 otherwise
    */
  int surfacePoint(const std::vector<double>& sKnots,
      const int sDegree,
      const std::vector<double>& tKnots,
      const int tDegree,
      const int Dim,
      const std::vector<double>& ctrlPts,
      const double s,
      const double t,
      std::vector<double>& outPts);

  //! Find a point on a curve
  /*!
    Uses deCasteljau to find the value of a Bezier curve
    at the given parametric value.
    \param ctrlPts The homogenous control points for this Bezier curve.
Note: polynomial degree is assumed from number of points.
\param Dim The dimension of the control points).
\param u Parametric value at which to evaluate the curve.  Should be
in [0,1].
\param pt Returns the coordinates of the point at u.
\return 0 if no errors, 1 otherwise.
*/
  int deCasteljauCurve(const std::vector<double>& ctrlPts,
      const int Dim,
      const double u,
      std::vector<double>& pt);


  //! Find a point on a surface
  /*!
    Uses deCasteljau to find the value of a Bezier surface
    at the given parametric values.
    \param ctrlPts The homogeneous control points for this Bezier surface.
    \param sDegree Polynomial degree in S.
    \param tDegree Polynomial degree in T.
    \param Dim The dimension of the control points.
    \param u Parametric value in s at which to evaluate the surface.  Should be
    in [0,1].
    \param v Parametric value in t at which to evaluate the surface.  Should be
    in [0,1].
    \param pt Returns the coordinates of the point at (u,v).
    \return 0 if no errors, 1 otherwise.
    */
  int deCasteljauSurface(const std::vector<double>& ctrlPts,
      const int sDegree,
      const int tDegree,
      const int Dim,
      const double u,
      const double v,
      std::vector<double>& pt);

  //! Find several points on a surface.
  /*!
    Uses deCasteljau to find the value of a Bezier surface
    at the given parametric values.
    \param ctrlPts The homogeneous control points for this Bezier surface.
    \param sDegree Polynomial degree in S.
    \param tDegree Polynomial degree in T.
    \param Dim The dimension of the curve control point.
    \param sVec Parametric values in s at which to evaluate the surface.  Should be
    in [0,1].
    \param tVec Parametric values in t at which to evaluate the surface.  Should be
    in [0,1].
    \param pts Returns the coordinates of the points at each (s,t) pair.
    \return 0 if no errors, 1 otherwise.
    */    
  int deCasteljauSurface(const std::vector<double>& ctrlPts,
      const int sDegree,
      const int tDegree,
      const int Dim,
      const std::vector<double>& sVec,
      const std::vector<double>& tVec,
      std::vector<double>& pts);

  //! Find points in a trivariate Bezier element
  /*!
    Uses deCasteljau to find the value of in a Bezier element
    at the given parametric values.
    \param ctrlPts The homogenous control points for this Bezier surface.
    \param sDegree Polynomial degree in S.
    \param tDegree Polynomial degree in T.
    \param uDegree Palynomial degree in U.
    \param sVec A list of parametric values in s at which to evaluate the surface.
    Should all be in [0,1].
    \param tVec A list of parametric values in t at which to evaluate the surface.  Should be
    in [0,1].
    \param uVec A list of parametric values in u at which to evaluate the surface.
    Should all be in [0,1].
    \param Dim The dimension of the control points.
    \param pts Returns the coordinates of the points at each (s,t,u) piont.  Points are
    ordered as usual in relation to s,t, and u.
    \return 0 if no errors, 1 otherwise.
    */
  int deCasteljauVolume(const std::vector<double>& ctrlPts,
      const int sDegree,
      const int tDegree,
      const int uDegree,
      const std::vector<double>& sVec,
      const std::vector<double>& tVec,
      const std::vector<double>& uVec,
      const int Dim,
      std::vector<double>& pts);

  //! Insert new knots into a curve
  /*!
    Computes a new curve from knot refinement.  Adapted from algorithm
    Piegl and Tiller "The NURBS Book" pg
    \param oldKnotVec The existing knot vector
    \param cDegree The polynomial degree of the curve
    \param Dim The dimension of the control points.
    \param oldCtrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param newKnots A vector of new knot values to insert into the old knot vector.
    \param newKnotVec the new knot vector.
    \param newCtrlPts the homogeneous control points corresponding to the new
    knot vector.
    \return 0 if no errors, 1 otherwise.
    */
  int knotRefinementCurve(const std::vector<double>& oldKnotVec,
      const int cDegree,
      const int Dim,
      const std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots,
      std::vector<double>& newKnotVec,
      std::vector<double>& newCtrlPts);

  //! Same as above except returns new vectors in place of old.    
  int knotRefinementCurve(std::vector<double>& oldKnotVec,
      const int cDegree,
      const int Dim,
      std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots);


  //! Insert new knots into a surface
  /*!
    Computes a new surface from knot refinement.  Adapted from
    algorithm in "The NURBS Book" pg 167.
    \param oldSKnotVec The existing knot vector in the S direction
    \param sDegree The polynomial degree in the S direction
    \param oldTKnotVec The existing knot vector in the T direction
    \param tDegree The polynomial degree in the T direction
    \param Dim The dimension of the control points.
    \param oldCtrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param newKnots A vector of new knot values to insert into the old knot vector.
    \param dir The knot vector to refine.  Should be 's' or 't' (case sensitive).
    \param newSKnotVec The new knot vector in the S direction.
    \param newTKnotVec The new knot vector in the T direction.
    \param newCtrlPts the homogeneous control points corresponding to the new
    knot vector.       
    \return 0 if no errors, 1 otherwise.
    */
  int knotRefinementSurface(const std::vector<double>& oldSKnotVec,
      const int sDegree,
      const std::vector<double>& oldTKnotVec,
      const int tDegree,
      const int Dim,
      const std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots,
      const char dir,
      std::vector<double>& newSKnotVec,
      std::vector<double>& newTKnotVec,
      std::vector<double>& newCtrlPts);

  //! Same as above except returns new vectors in place of old.
  int knotRefinementSurface(std::vector<double>& oldSKnotVec,
      const int sDegree,
      std::vector<double>& oldTKnotVec,
      const int tDegree,
      const int Dim,
      std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots,
      const char dir);

  //! Insert new knots into a volume
  /*!
    Computes a new volume from knot refinement.
    \param oldSKnotVec The existing knot vector in the S direction
    \param sDegree The polynomial degree in the S direction
    \param oldTKnotVec The existing knot vector in the T direction
    \param tDegree The polynomial degree in the T direction
    \param oldUKnotVec the existing knot vector in the U direction
    \param uDegree The polynomial degree in the U direction
    \param oldCtrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param newKnots A vector of new knot values to insert into the old knot vector.
    \param dir The knot vector to refine.  Should be 's', 't', or, 'u'  (case sensitive).
    \param Dim The dimension of the control points.
    \param newSKnotVec The new knot vector in the S direction.
    \param newTKnotVec The new knot vector in the T direction.
    \param newUKnotVec the new knot vector in the U direction.
    \param newCtrlPts the homogeneous control points corresponding to the new
    knot vector.       
    \return 0 if no errors, 1 otherwise.
    */
  int knotRefinementVolume(const std::vector<double>& oldSKnotVec,
      const int sDegree,
      const std::vector<double>& oldTKnotVec,
      const int tDegree,
      const std::vector<double>& oldUKnotVec,
      const int uDegree,
      const std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots,
      const char dir,
      const int Dim,
      std::vector<double>& newSKnotVec,
      std::vector<double>& newTKnotVec,
      std::vector<double>& newUKnotVec,
      std::vector<double>& newCtrlPts);

  //! Same as above except returns new vectors in place of old.
  int knotRefinementVolume(std::vector<double>& oldSKnotVec,
      const int sDegree,
      std::vector<double>& oldTKnotVec,
      const int tDegree,
      std::vector<double>& oldUKnotVec,
      const int uDegree,
      std::vector<double>& oldCtrlPts,
      const std::vector<double>& newKnots,
      const char dir,
      const int Dim);

  //! Degree elevate a curve
  /*!
    Computes a new curve from knot refinement.  Adapted from algorithm
    Piegl and Tiller "The NURBS Book" pg. 206.
    \param oldKnotVec The existing knot vector
    \param cDegree The polynomial degree of the curve
    \param Dim The dimension of the control points.
    \param oldCtrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param addDegrees The number of times to degree elevate the curve.
    \param newKnotVec the new knot vector.
    \param newCtrlPts the homogeneous control points corresponding to the new
    knot vector.
    \return 0 if no errors, 1 otherwise.
    */    
  int degreeElevateCurve(const std::vector<double>& oldKnotVec,
      const int cDegree,
      const int Dim,
      const std::vector<double>& oldCtrlPts,
      const int addDegrees,
      std::vector<double>& newKnotVec,
      std::vector<double>& newCtrlPts);

  //! Same as above except overwrite old vectors with new values.    
  int  degreeElevateCurve(std::vector<double>& oldKnotVec,
      const int cDegree,
      const int Dim,
      std::vector<double>& oldCtrlPts,
      const int addDegrees);

  //! Degree elevate a surface
  /*!
    Computes a new curve from knot refinement.  Adapted from algorithm
    Piegl and Tiller "The NURBS Book" pg. 206.
    \param oldSKnotVec The existing knot vector in the S direction
    \param sDegree The polynomial degree in the S direction
    \param oldTKnotVec The existing knot vector in the T direction
    \param tDegree The polynomial degree in the T direction
    \param Dim The dimension of the control points.
    \param oldCtrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param addDegrees The number of times to degree elevate the curve.
    \param dir The knot vector to refine.  Should be 's' or 't' (case sensitive).
    \param newSKnotVec The new knot vector in the S direction.
    \param newTKnotVec The new knot vector in the T direction.
    \param newCtrlPts the homogeneous control points corresponding to the new
    knot vector.
    \return 0 if no errors, 1 otherwise.
    */    
  int degreeElevateSurface(const std::vector<double>& oldSKnotVec,
      const int sDegree,
      const std::vector<double>& oldTKnotVec,
      const int tDegree,
      const int Dim,
      const std::vector<double>& oldCtrlPts,
      const int addDegrees,
      const char dir,
      std::vector<double>& newSKnotVec,
      std::vector<double>& newTKnotVec,
      std::vector<double>& newCtrlPts);

  //! Same as above except overwrite old vectors with new values.
  int degreeElevateSurface(std::vector<double>& oldSKnotVec,
      const int sDegree,
      std::vector<double>& oldTKnotVec,
      const int tDegree,
      const int Dim,
      std::vector<double>& oldCtrlPts,
      const int addDegrees,
      const char dir);

  //! Degree elevate a volume 
  /*! 
    Computes a new curve from knot refinement.  Adapted from algorithm 
    Piegl and Tiller "The NURBS Book" pg. 206. 
    \param oldSKnotVec The existing knot vector in the S direction 
    \param sDegree The polynomial degree in the S direction 
    \param oldTKnotVec The existing knot vector in the T direction 
    \param tDegree The polynomial degree in the T direction 
    \param oldUKnotVec The existing knot vector in the U direction 
    \param uDegree The polynomial degree in the U direction 
    \param Dim The dimension of the control points. 
    \param oldCtrlPts The homogeneous control points corresponding to the old 
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}. 
    \param addDegrees The number of times to degree elevate the curve. 
    \param dir The knot vector to refine.  Should be 's', 't', or 'u' (case sensitive). 
    \param newSKnotVec The new knot vector in the S direction. 
    \param newTKnotVec The new knot vector in the T direction. 
    \param newUKnotVec The new knot vector in the U direction. 
    \param newCtrlPts the homogeneous control points corresponding to the new 
    knot vector. 
    \return 0 if no errors, 1 otherwise. 
    */     
  int degreeElevateVolume(const std::vector<double>& oldSKnotVec, 
      const int sDegree, 
      const std::vector<double>& oldTKnotVec, 
      const int tDegree, 
      const std::vector<double>& oldUKnotVec, 
      const int uDegree, 
      const int Dim, 
      const std::vector<double>& oldCtrlPts, 
      const int addDegrees, 
      const char dir, 
      std::vector<double>& newSKnotVec, 
      std::vector<double>& newTKnotVec, 
      std::vector<double>& newUKnotVec, 
      std::vector<double>& newCtrlPts); 

  //! Same as above except overwrite old vectors with new values. 
  int degreeElevateVolume(std::vector<double>& oldSKnotVec, 
      const int sDegree, 
      std::vector<double>& oldTKnotVec, 
      const int tDegree, 
      std::vector<double>& oldUKnotVec, 
      const int uDegree, 
      const int Dim, 
      std::vector<double>& oldCtrlPts, 
      const int addDegrees, 
      const char dir); 

  //! Decomposes a curve into Bezier segments
  /*!
    Adapted from algorithm in "The NURBS Book" by Piegl and Tilelr pg 173.
    \param knotVec The knot vector for the curve
    \param cDegree The polynomial degree of the curve
    \param ctrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.
    \param Dim The dimension of the control points.
    \param bCtrlPts The homogeneous control points for the Bezier segments.
    The points are stored the same as for ctrlPts.
    \param bSegments A list of the new Bezier Segments.  Calling function
    is responsible for deleting the segments.
    \param calcbCtrlPts Set to true to calculate the control points for the
    Bezier segments.  If set to false bCtrlPts will be empty at return
    and ctrlPts are ignored.       
    \return 0 if no errors, 1 otherwise.
    */
  int decomposeCurve(const std::vector<double>& knotVec,
      const int cDegree,
      const std::vector<double>& ctrlPts,
      const int Dim,
      std::vector<double>& bCtrlPts,
      std::vector<BezierElem*>& bSegments,
      const bool calcbCtrlPts = true);

  //! Decomposes a surface into Bezier segments
  /*!
    Adapted from algorithm in "The NURBS Book" by Piegl and Tilelr pg 177.
    \param sKnots The knot vector for the S direction
    \param sDegree The polynomial degree in the S direction
    \param tKnots The knot vector in the T direction
    \param tDegree the polynomial degree in the T direction
    \param ctrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.       
    \param Dim The dimension of the control points.
    \param bCtrlPts The homogeneous control points for the Bezier segments.
    The points are stored the same as for ctrlPts.
    \param bSegments A list of the new Bezier Segments.  Calling function
    is responsible for deleting the segments.
    \param buildFullOps If set to true then the Bezier elements are built with
    full extraction operators initialized.  If set to false then
    the operators are left in component form.
    \return 0 if no errors, 1 otherwise.
    */
  int decomposeSurface(const std::vector<double>& sKnots,
      const int sDegree,
      const std::vector<double>& tKnots,
      const int tDegree,
      const std::vector<double>& ctrlPts,
      const int Dim,
      std::vector<double>& bCtrlPts,
      std::vector<BezierElem*>& bSegments,
      const bool buildFullOps = true);

  //! Decomposes a surface into Bezier segments
  /*!
    Adapted from algorithm in "The NURBS Book" by Piegl and Tilelr pg 177.
    \param sKnots The knot vector for the S direction
    \param sDegree The polynomial degree in the S direction
    \param tKnots The knot vector in the T direction
    \param tDegree the polynomial degree in the T direction
    \param uKnots The knot vector for the U direction
    \param uDegree The polynomial degree int the U direction
    \param ctrlPts The homogeneous control points corresponding to the old
    knot vector.  Points stored as {x1,y1,z1,w1,...,xn,yn,zn,wn}.       
    \param bCtrlPts The homogeneous control points for the Bezier segments.
    The points are stored the same as for ctrlPts.
    \param Dim The dimension of the control points.      
    \param bSegments A list of the new Bezier Segments.  Calling function
    is responsible for deleting the segments.
    \param buildFullOps If set to true (default) then the Bezier elements are built with
    full extraction operators initialized.  If set to false then
    the operators are left in component form.
    \return 0 if no errors, 1 otherwise.
    */
  int decomposeVolume(const std::vector<double>& sKnots,
      const int sDegree,
      const std::vector<double>& tKnots,
      const int tDegree,
      const std::vector<double>& uKnots,
      const int uDegree,
      const std::vector<double>& ctrlPts,
      const int Dim,
      std::vector<double>& bCtrlPts,
      std::vector<BezierElem*>& bSegments,
      const bool buildFullOps = true);

  //! Use the Coefficient matrices of the Bezier segments to tranform points.
  /*!
    This function transforms global control points to control points for the
    Bezier segments using the coefficient matrices.
    \param inPts The homogeneous coordinates of the points that will be transformed.
Note: The number and order (but not dimension) of the points must match
that of the points used to create the bezier decomposition.
\param bSegments The list of Bezier segments
\param Dim the dimension of the control points.
\param outPts The homogeneous coordinates of the resulting points.
\return 0 if no errors, 1 otherwise.
*/
  int transformPoints(const std::vector<double>& inPts,
      const std::vector<BezierElem*>& bSegments,
      const int Dim,
      std::vector<double>& outPts);

  /////////////////////////////////////////////////////////////////
  ////  Place Helper Functions below this line ////
  /////////////////////////////////////////////////////////////////
  int binomialCoefficient(const int n, const int m);


}
#endif

//EOF

