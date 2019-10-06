#include <math.h>
#include <iostream>
#include <assert.h>
#include <map>

#include "NURBS_Tools.hpp"
#include "NURBS_Bezier.hpp"

int NURBS_T::projectUp(std::vector<double>& ctrlPts,
                       const int spatialDim)
{
  const int cpDim = spatialDim+1;
  const int numPts = ctrlPts.size()/cpDim;
  
  int ii, jj;
  for( ii = 0; ii < numPts; ii++ )
  {
    double w = ctrlPts[cpDim*ii+spatialDim];
    for( jj = 0; jj < spatialDim; jj++ )
        ctrlPts[cpDim*ii+jj] = ctrlPts[cpDim*ii+jj] * w;
  }
    
  return 0;
}

int NURBS_T::projectDown(std::vector<double>& ctrlPts,
                         const int spatialDim)
{
  const int cpDim = spatialDim+1;
  const int numPts = ctrlPts.size()/cpDim;
  
  int ii, jj;
  for( ii = 0; ii < numPts; ii++ )
  {
    double w = ctrlPts[cpDim*ii+spatialDim];
    for( jj = 0; jj < spatialDim; jj++ )
        ctrlPts[cpDim*ii+jj] = ctrlPts[cpDim*ii+jj] / w;
  }
  
  return 0;
}

int NURBS_T::findSpan(const std::vector<double>& knotVec,
                      const int cDegree,
                      const double value)
{
  if( value < knotVec[0] || value > knotVec.back() )
  {
    std::cerr << "Error: value outside knot vector end points in findSpan()" << std::endl;
    return -1;
  }

  int high = (knotVec.size()-1) - cDegree;
  int low = cDegree;
  int mid = (low + high)/2;
  
    //Special end case.
  if( value == knotVec[high] )
      return high - 1;
  
  while( value < knotVec[mid] || value >= knotVec[mid+1] )
  {
    if( value < knotVec[mid] )
        high = mid;
    else
        low = mid;
    mid = (high + low)/2;
  }
  return mid;
}

int NURBS_T::knotRefinementCurve(const std::vector<double>& oldKnotVec,
                                 const int cDegree,
                                 const int Dim,
                                 const std::vector<double>& oldCtrlPts,
                                 const std::vector<double>& newKnots,
                                 std::vector<double>& newKnotVec,
                                 std::vector<double>& newCtrlPts)
{
    //Make sure the new lists are clean.
  newKnotVec.clear();
  newCtrlPts.clear();

  int beginSpan = findSpan(oldKnotVec, cDegree, newKnots.front());
  int endSpan = findSpan(oldKnotVec, cDegree, newKnots.back() );

  if( beginSpan < 0 || endSpan < 0 )
      return 1;

  endSpan++;
  
    //Save unaltered control points.
  int ii, jj;
  const int  cpDim = Dim;
  newCtrlPts.resize( oldCtrlPts.size() + cpDim*newKnots.size(), 0.0 );
  for( ii = 0; ii <= beginSpan - cDegree; ii++ )
  {
    for( jj = 0; jj < cpDim; jj++ )
        newCtrlPts[cpDim*ii+jj] = oldCtrlPts[cpDim*ii+jj];
  }
  for( ii = endSpan-1; ii < (int)oldCtrlPts.size()/cpDim; ii++ )
  {
    for( jj = 0; jj < cpDim; jj++ )
        newCtrlPts[cpDim*(ii+newKnots.size())+jj] = oldCtrlPts[cpDim*ii+jj];
  }

    //Load the new knot vector with knots that will not change.
  newKnotVec.resize(oldKnotVec.size()+newKnots.size(), 0.0 );
  for( ii = 0; ii <= beginSpan; ii++ )
      newKnotVec[ii] = oldKnotVec[ii];
  for( ii = endSpan+cDegree; ii < (int)oldKnotVec.size(); ii++ )
      newKnotVec[ii+newKnots.size()] = oldKnotVec[ii];

  
  int kk,mm;
  ii = endSpan+cDegree-1;
  kk = endSpan+cDegree+newKnots.size()-1;

  for( jj = newKnots.size()-1; jj >= 0; jj-- )
  {
    while( newKnots[jj] <= oldKnotVec[ii] && ii > beginSpan )
    {
      for( mm = 0; mm < cpDim; mm++ )
          newCtrlPts[cpDim*(kk-cDegree-1)+mm] = oldCtrlPts[cpDim*(ii-cDegree-1)+mm];

      newKnotVec[kk] = oldKnotVec[ii];
      kk--;
      ii--;
    }
    for( mm = 0; mm < cpDim; mm++ )
          newCtrlPts[cpDim*(kk-cDegree-1)+mm] = newCtrlPts[cpDim*(kk-cDegree)+mm];

    int ll;
    for( ll = 1; ll <= cDegree; ll++ )
    {
      int ind = kk-cDegree+ll;
      double alpha = newKnotVec[kk+ll] - newKnots[jj];
      if( fabs(alpha) == 0.0 )
      {
        for( mm = 0; mm < cpDim; mm++ )
          newCtrlPts[cpDim*(ind-1)+mm] = newCtrlPts[cpDim*(ind)+mm];
      }
      else
      {
        alpha = alpha / (newKnotVec[kk+ll]-oldKnotVec[ii-cDegree+ll]);
        for( mm = 0; mm < cpDim; mm++ )
            newCtrlPts[cpDim*(ind-1)+mm] =  alpha*newCtrlPts[cpDim*(ind-1)+mm]
                + (1.0-alpha)*newCtrlPts[cpDim*ind+mm];
      }
    }
    newKnotVec[kk] = newKnots[jj];
    kk--;
  }
  
  return 0;
}

int NURBS_T::knotRefinementCurve(std::vector<double>& oldKnotVec,
                                 const int cDegree,
                                 const int Dim,
                                 std::vector<double>& oldCtrlPts,
                                 const std::vector<double>& newKnots)
{
  std::vector<double> newKnotVec;
  std::vector<double> newPts;

  if( 0 != knotRefinementCurve( oldKnotVec, cDegree, Dim,
                                oldCtrlPts, newKnots,
                                newKnotVec, newPts ) )
      return 1;

  oldKnotVec.clear();
  oldCtrlPts.clear();
  oldKnotVec.insert( oldKnotVec.end(), newKnotVec.begin(), newKnotVec.end() );
  oldCtrlPts.insert( oldCtrlPts.end(), newPts.begin(), newPts.end() );
  
  return 0;
}

int NURBS_T::knotRefinementSurface(const std::vector<double>& oldUKnotVec,
                                   const int uDegree,
                                   const std::vector<double>& oldVKnotVec,
                                   const int vDegree,
                                   const int Dim,
                                   const std::vector<double>& oldCtrlPts,
                                   const std::vector<double>& newKnots,
                                   const char dir,
                                   std::vector<double>& newUKnotVec,
                                   std::vector<double>& newVKnotVec,
                                   std::vector<double>& newCtrlPts)
{
  newUKnotVec.clear();
  newVKnotVec.clear();
  newCtrlPts.clear();
  const int  cpDim = Dim;
  std::vector<double>::iterator inIter;
  std::vector<double>::const_iterator outIter;
  
  if( 's' == dir )
  {
    int beginSpan = findSpan(oldUKnotVec, uDegree, newKnots.front());
    int endSpan = findSpan(oldUKnotVec, uDegree, newKnots.back() );

    if( beginSpan < 0 || endSpan < 0 )
        return 1;
    
    endSpan++;

    int ii, jj;
    
      //Load the new knot vectors with knots that will not change.
    newUKnotVec.resize(oldUKnotVec.size()+newKnots.size(), 0.0 );
    for( ii = 0; ii <= beginSpan; ii++ )
        newUKnotVec[ii] = oldUKnotVec[ii];
    for( ii = endSpan+uDegree; ii < (int)oldUKnotVec.size(); ii++ )
        newUKnotVec[ii+newKnots.size()] = oldUKnotVec[ii];

    newVKnotVec.insert( newVKnotVec.end(), oldVKnotVec.begin(), oldVKnotVec.end() );

      //Save unaltered control points.
    int row;
    const int uCols = oldUKnotVec.size() - uDegree - 1;
    const int vRows = oldVKnotVec.size() - vDegree - 1;
    const int newUCols = uCols+newKnots.size();
    newCtrlPts.resize( oldCtrlPts.size()+(vRows*cpDim*newKnots.size()), 0.0 );
    
    for( row = 0; row < vRows; row++ )
    {
      for( ii = 0; ii <= beginSpan - uDegree; ii++ )
      {
        for( jj = 0; jj < cpDim; jj++ )
            newCtrlPts[cpDim*(row*newUCols+ii)+jj]
                = oldCtrlPts[cpDim*(row*uCols+ii)+jj];
      }
      for( ii = endSpan-1; ii < uCols; ii++ )
      {
        for( jj = 0; jj < cpDim; jj++ )
            newCtrlPts[cpDim*(row*newUCols+ii+newKnots.size())+jj]
                = oldCtrlPts[cpDim*(row*uCols+ii)+jj];
      }
    }
    
      //Calculate new control points.
    int kk,mm;
    ii = endSpan+uDegree-1;
    kk = endSpan+uDegree+newKnots.size()-1;

    for( jj = newKnots.size()-1; jj >= 0; jj-- )
    {
      while( newKnots[jj] <= oldUKnotVec[ii] && ii > beginSpan )
      {
        newUKnotVec[kk] = oldUKnotVec[ii];
        for( row = 0; row < vRows; row++ )
        {
          for( mm = 0; mm < cpDim; mm++ )
          {
            newCtrlPts[cpDim*(row*newUCols+kk-uDegree-1)+mm]
                = oldCtrlPts[cpDim*(row*uCols+ii-uDegree-1)+mm];
          }
        }
        kk--;
        ii--;
      }
      for( row = 0; row < vRows; row++ )
      {
        for( mm = 0; mm < cpDim; mm++ )
        {
          newCtrlPts[cpDim*(row*newUCols+kk-uDegree-1)+mm]
              = newCtrlPts[cpDim*(row*newUCols+kk-uDegree)+mm];
        }
      }

      int ll;
      for( ll = 1; ll <= uDegree; ll++ )
      {
        int ind = kk-uDegree+ll;
        double alpha = newUKnotVec[kk+ll] - newKnots[jj];
        if( fabs(alpha) == 0.0 )
        {
          for( row = 0; row < vRows; row++ )
          {
            for( mm = 0; mm < cpDim; mm++ )
            {
              newCtrlPts[cpDim*(row*newUCols+ind-1)+mm]
                  = newCtrlPts[cpDim*(row*newUCols+ind)+mm];
            }
          }
        }
        else
        {
          alpha = alpha / (newUKnotVec[kk+ll]-oldUKnotVec[ii-uDegree+ll]);
          for( row = 0; row < vRows; row++ )
          {
            for( mm = 0; mm < cpDim; mm++ )
            {
              newCtrlPts[cpDim*(row*newUCols+ind-1)+mm] =
                  alpha*newCtrlPts[cpDim*(row*newUCols+ind-1)+mm]
                  + (1.0-alpha)*newCtrlPts[cpDim*(row*newUCols+ind)+mm];
            }
          }
        }
      }
      newUKnotVec[kk] = newKnots[jj];
      kk--;
    }
  }
  else if( 't' == dir )
  {
    int beginSpan = findSpan(oldVKnotVec, vDegree, newKnots.front());
    int endSpan = findSpan(oldVKnotVec, vDegree, newKnots.back() );

    if( beginSpan < 0 || endSpan < 0 )
        return 1;
    
    endSpan++;

    int ii, jj;
    
      //Load the new knot vectors with knots that will not change.
    newVKnotVec.resize(oldVKnotVec.size()+newKnots.size(), 0.0 );
    for( ii = 0; ii <= beginSpan; ii++ )
        newVKnotVec[ii] = oldVKnotVec[ii];
    for( ii = endSpan+vDegree; ii < (int)oldVKnotVec.size(); ii++ )
        newVKnotVec[ii+newKnots.size()] = oldVKnotVec[ii];

    newUKnotVec.insert( newUKnotVec.end(), oldUKnotVec.begin(), oldUKnotVec.end() );

      //Save unaltered control points.
    int col;
    const int uCols = oldUKnotVec.size() - uDegree - 1;
    const int vRows = oldVKnotVec.size() - vDegree - 1;
    newCtrlPts.resize( oldCtrlPts.size()+(uCols*cpDim*newKnots.size()), 0.0 );
    
    for( col = 0; col < uCols; col++ )
    {
      for( ii = 0; ii <= beginSpan - vDegree; ii++ )
      {
        for( jj = 0; jj < cpDim; jj++ )
            newCtrlPts[cpDim*(col+uCols*ii)+jj]
                = oldCtrlPts[cpDim*(col+uCols*ii)+jj];
      }
      for( ii = endSpan-1; ii < vRows; ii++ )
      {
        for( jj = 0; jj < cpDim; jj++ )
            newCtrlPts[cpDim*(col+uCols*(ii+newKnots.size()))+jj]
                = oldCtrlPts[cpDim*(col+uCols*ii)+jj];
      }
    }
    
      //Calculate new control points.
    int kk,mm;
    ii = endSpan+vDegree-1;
    kk = endSpan+vDegree+newKnots.size()-1;

    for( jj = newKnots.size()-1; jj >= 0; jj-- )
    {
      while( newKnots[jj] <= oldVKnotVec[ii] && ii > beginSpan )
      {
        newVKnotVec[kk] = oldVKnotVec[ii];
        for( col = 0; col < uCols; col++ )
        {
          for( mm = 0; mm < cpDim; mm++ )
          {
            newCtrlPts[cpDim*(col+uCols*(kk-vDegree-1))+mm]
                = oldCtrlPts[cpDim*(col+uCols*(ii-vDegree-1))+mm];
          }
        }
        kk--;
        ii--;
      }
      for( col = 0; col < uCols; col++ )
      {
        for( mm = 0; mm < cpDim; mm++ )
        {
          newCtrlPts[cpDim*(col+uCols*(kk-vDegree-1))+mm]
              = newCtrlPts[cpDim*(col+uCols*(kk-vDegree))+mm];
        }
      }

      int ll;
      for( ll = 1; ll <= vDegree; ll++ )
      {
        int ind = kk-vDegree+ll;
        double alpha = newVKnotVec[kk+ll] - newKnots[jj];
        if( fabs(alpha) == 0.0 )
        {
          for( col = 0; col < uCols; col++ )
          {
            for( mm = 0; mm < cpDim; mm++ )
            {
              newCtrlPts[cpDim*(col+uCols*(ind-1))+mm]
                  = newCtrlPts[cpDim*(col+uCols*ind)+mm];
            }
          }
        }
        else
        {
          alpha = alpha / (newVKnotVec[kk+ll]-oldVKnotVec[ii-vDegree+ll]);
          for( col = 0; col < uCols; col++ )
          {
            for( mm = 0; mm < cpDim; mm++ )
            {
              newCtrlPts[cpDim*(col+uCols*(ind-1))+mm] =
                  alpha*newCtrlPts[cpDim*(col+uCols*(ind-1))+mm]
                  + (1.0-alpha)*newCtrlPts[cpDim*(col+uCols*ind)+mm];
            }
          }
        }
      }
      newVKnotVec[kk] = newKnots[jj];
      kk--;
    }    
  }
  else
  {
    std::cerr << "Incorrect direction parameter given in NURBS_T::knotRefinementSurface()"
              << std::endl;
    return 1;
  }
  
  return 0;
}

int NURBS_T::knotRefinementSurface(std::vector<double>& oldUKnotVec,
                                   const int uDegree,
                                   std::vector<double>& oldVKnotVec,
                                   const int vDegree,
                                   const int Dim,
                                   std::vector<double>& oldCtrlPts,
                                   const std::vector<double>& newKnots,
                                   const char dir)
{
  std::vector<double> newUKnots;
  std::vector<double> newVKnots;
  std::vector<double> newPts;

  if( 0 != knotRefinementSurface( oldUKnotVec, uDegree, oldVKnotVec, vDegree,
                                  Dim, oldCtrlPts, newKnots, dir,
                                  newUKnots, newVKnots, newPts ) )
      return 1;

  oldUKnotVec.clear();
  oldVKnotVec.clear();
  oldCtrlPts.clear();
  oldUKnotVec.insert(oldUKnotVec.end(), newUKnots.begin(), newUKnots.end() );
  oldVKnotVec.insert(oldVKnotVec.end(), newVKnots.begin(), newVKnots.end() );
  oldCtrlPts.insert(oldCtrlPts.end(), newPts.begin(), newPts.end() );

  return 0;
}

int NURBS_T::knotRefinementVolume(const std::vector<double>& oldSKnotVec,
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
                                  std::vector<double>& newCtrlPts)
{
  const int cpDim = Dim;
  
  newSKnotVec.clear();
  newTKnotVec.clear();
  newUKnotVec.clear();
  newCtrlPts.clear();

  const int sCols = oldSKnotVec.size() - sDegree - 1;
  const int tRows = oldTKnotVec.size() - tDegree - 1;
  const int uLayers = oldUKnotVec.size() - uDegree - 1;
  std::vector<double> tempInPts;
  std::vector<double> tempOutPts;
    
    //Do surface knot refinement on each layer.
  if( 's' == dir || 't' == dir )
  {
      //Save unchanged u knots.
    newUKnotVec.insert(newUKnotVec.end(), oldUKnotVec.begin(), oldUKnotVec.end() );
    
    int layer;
    for( layer = 0; layer < uLayers; layer++ )
    {
        //Get the control points for the current layer.
      int row, col;
      for( row = 0; row < tRows; row++ )
      {
        for( col = 0; col < sCols; col++ )
        {
          int ii;
          for( ii = 0; ii < Dim; ii++ )
              tempInPts.push_back( oldCtrlPts[Dim*(layer*sCols*tRows+row*sCols+col)+ii] );
        }
      }
      
        //Call knot refinement on the "surface".
      if( 0 != knotRefinementSurface(oldSKnotVec, sDegree, oldTKnotVec, tDegree,
                                     cpDim, tempInPts, newKnots, dir,
                                     newSKnotVec, newTKnotVec, tempOutPts) )
          return 1;

      newCtrlPts.insert( newCtrlPts.end(), tempOutPts.begin(), tempOutPts.end() );

      tempInPts.clear();
      tempOutPts.clear();
    }
  }
  else if( 'u' == dir )
  {
      //Save unchanged t knots.
    newTKnotVec.insert(newTKnotVec.end(), oldTKnotVec.begin(), oldTKnotVec.end() );

      //Resize the new control point vector.
    newCtrlPts.resize( oldCtrlPts.size() + (Dim*sCols*tRows*newKnots.size()), 0.0 );

    const int newULayers = uLayers + newKnots.size();
    
    int col, row, layer;
    for( row = 0; row < tRows; row++ )
    {
        //Get the control points for the current layer.      
      for( layer = 0; layer < uLayers; layer++ )
      {
        for( col = 0; col < sCols; col++ )
        {
          int ii;
          for( ii = 0; ii < Dim; ii++ )
              tempInPts.push_back( oldCtrlPts[Dim*(layer*sCols*tRows+row*sCols+col)+ii] );
        }
      }

        //Call knot refinement on the "surface".
      if( 0 != knotRefinementSurface(oldSKnotVec, sDegree, oldUKnotVec, uDegree,
                                     cpDim, tempInPts, newKnots, 't',
                                     newSKnotVec, newUKnotVec, tempOutPts) )
          return 1;

      std::vector<double>::iterator iter = tempOutPts.begin();
      for( layer = 0; layer < newULayers; layer++ )
      {
        for( col = 0; col < sCols; col++ )
        {
          int ii;
          for( ii = 0; ii < Dim; ii++ )
          {
            newCtrlPts[Dim*(layer*sCols*tRows+row*sCols+col)+ii] = *iter;
            iter++;
          }
        }
      }

      tempInPts.clear();
      tempOutPts.clear();
    }
  }
  else
  {
    std::cerr << "Incorrect direction parameter given in NURBS_T::knotRefinementVolume()"
              << std::endl;
    return 1;
  }
  
  return 0;
}

int NURBS_T::knotRefinementVolume(std::vector<double>& oldSKnotVec,
                                  const int sDegree,
                                  std::vector<double>& oldTKnotVec,
                                  const int tDegree,
                                  std::vector<double>& oldUKnotVec,
                                  const int uDegree,
                                  std::vector<double>& oldCtrlPts,
                                  const std::vector<double>& newKnots,
                                  const char dir,
                                  const int Dim)
{
  std::vector<double> newSKnots;
  std::vector<double> newTKnots;
  std::vector<double> newUKnots;
  std::vector<double> newPts;

  if( 0 != knotRefinementVolume( oldSKnotVec, sDegree,
                                 oldTKnotVec, tDegree,
                                 oldUKnotVec, uDegree,
                                 oldCtrlPts, newKnots, dir, Dim,
                                 newSKnots, newTKnots, newUKnots, newPts ) )
      return 1;

  oldSKnotVec.clear();
  oldTKnotVec.clear();
  oldUKnotVec.clear();
  oldCtrlPts.clear();
  oldSKnotVec.insert(oldSKnotVec.end(), newSKnots.begin(), newSKnots.end() );
  oldTKnotVec.insert(oldTKnotVec.end(), newTKnots.begin(), newTKnots.end() );
  oldUKnotVec.insert(oldUKnotVec.end(), newUKnots.begin(), newUKnots.end() );
  oldCtrlPts.insert(oldCtrlPts.end(), newPts.begin(), newPts.end() );
  
  return 0;
}

int NURBS_T::deCasteljauCurve(const std::vector<double>& ctrlPts,
                              const int Dim,
                              const double u,
                              std::vector<double>& pt)
{
    //Use local copy of points
  std::vector<double> pts;
  pts.insert( pts.end(), ctrlPts.begin(), ctrlPts.end() );

  const int cpDim = Dim;
  const int numPts = ctrlPts.size() / cpDim;
  
  int ii, kk, mm;
  for( kk = 1; kk < numPts; kk++ )
  {
    for( ii = 0; ii < numPts-kk; ii++ )
    {
      for( mm = 0; mm < cpDim; mm++ )
          pts[cpDim*ii+mm]=(1.0-u)*pts[cpDim*ii+mm] + u*pts[cpDim*(ii+1)+mm];
    }
  }

  pt.clear();
  for( ii = 0; ii < cpDim; ii++ )
      pt.push_back( pts[ii] );
  
  return 0;
}

int NURBS_T::deCasteljauSurface(const std::vector<double>& ctrlPts,
                                const int sDegree,
                                const int tDegree,
                                const int Dim,
                                const double u,
                                const double v,
                                std::vector<double>& pt)
{
  pt.clear();
  
  std::vector<double> pts;

  const int cpDim = Dim;
  if( sDegree <= tDegree )
  {
      //Find column points first.
    std::vector<double> locPts;
    int ii;
    std::vector<double>::const_iterator iter = ctrlPts.begin();
    for( ii = 0; ii <= tDegree; ii++ )
    {
        //Get the ii row of points.
      pts.clear();
      pts.insert( pts.end(), iter, iter+cpDim*(tDegree+1) );
      iter += cpDim*(tDegree+1);

      if( 0 != deCasteljauCurve( pts, Dim, u, pt ) )
      {
        std::cerr << "ERROR: in NURBS_T::deCasteljauSurface()" << std::endl;
        return 1;
      }
      locPts.insert( locPts.end(), pt.begin(), pt.end() );

      pt.clear();
    }

    if( 0 != deCasteljauCurve( locPts, Dim, v, pt ) )
    {
      std::cerr << "ERROR: in NURBS_T::deCasteljauSurface() -- 1" << std::endl;
      return 1;
    }    
  }
  else
  {
      //Find row points first.
    std::vector<double> locPts;
    std::vector<double>::const_iterator iter;
    
    int ii;
    for( ii = 0; ii <= sDegree; ii++ )
    {
      pts.clear();

        //Get the ii column of points from the from the jj row. 
      int jj;
      for( jj = 0; jj <= tDegree; jj++ )
      {
        iter = ctrlPts.begin() + cpDim*(jj*(sDegree+1)+ii);
        int mm;
        for( mm = 0; mm < cpDim; mm++ )
        {
          pts.push_back( *iter );
          iter++;
        }
      }

      if( 0 != deCasteljauCurve( pts, Dim, v, pt ) )
      {
        std::cerr << "ERROR: in NURBS_T::deCasteljauSurface()" << std::endl;
        return 1;
      }
      locPts.insert( locPts.end(), pt.begin(), pt.end() );

      pt.clear();
    }

    if( 0 != deCasteljauCurve( locPts, Dim, u, pt ) )
    {
      std::cerr << "ERROR: in NURBS_T::deCasteljauSurface() -- 2" << std::endl;
      return 1;
    }        
    
  }
  
  return 0;
}

int NURBS_T::deCasteljauVolume(const std::vector<double>& ctrlPts,
                               const int sDegree,
                               const int tDegree,
                               const int uDegree,
                               const std::vector<double>& sVec,
                               const std::vector<double>& tVec,
                               const std::vector<double>& uVec,
                               const int Dim,
                               std::vector<double>& pts)
{
  pts.clear();

  std::vector<double> surfPts;
  std::vector<double> tempPts;
  std::vector<double> locPts;
  std::vector<double> pt;
  const int cpDim = Dim;
  if( sDegree <= tDegree && sDegree <= uDegree )
  {    
    pts.resize(cpDim*sVec.size()*tVec.size()*uVec.size(), 0);
    std::vector<double>::iterator p_iter;
    
    int ii;
    for( ii = 0; ii < (int)sVec.size(); ii++ )
    {
      double s = sVec[ii];

        //Get the t-u surface at this s value.
      int jj, kk;
      std::vector<double>::const_iterator iter = ctrlPts.begin();
      for( kk = 0; kk <= uDegree; kk++ )
      {
        for( jj = 0; jj <= tDegree; jj++ )
        {
          locPts.insert(locPts.end(), iter, iter+cpDim*(sDegree+1) );
          iter += cpDim*(sDegree+1);

          if( 0 != deCasteljauCurve( locPts, cpDim, s, pt ) )
          {
            std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 1" << std::endl;
            return 1;
          }
          tempPts.insert( tempPts.end(), pt.begin(), pt.end() );
          
          pt.clear();
          locPts.clear();
        }
      }

        //Now get the points from the surface.
      if( 0 != deCasteljauSurface( tempPts, tDegree, uDegree, cpDim, tVec, uVec, surfPts ) )
      {
        std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 2" << std::endl;
        return 1;
      }

      std::vector<double>::iterator s_iter = surfPts.begin();
      p_iter = pts.begin() + cpDim*ii;
      for( kk = 0; kk < (int)uVec.size(); kk++ )
      {
        for( jj = 0; jj < (int)tVec.size(); jj++ )
        {
          int mm;
          for( mm = 0; mm < cpDim; mm++ )
          {
            *(p_iter + mm) = *s_iter;
            s_iter++;
          }
          
          p_iter += cpDim*(sVec.size());
        }
      }
      
      tempPts.clear();
      surfPts.clear();
    }
  }
  if( tDegree <= sDegree && tDegree <= uDegree )
  {
    pts.resize(cpDim*sVec.size()*tVec.size()*uVec.size(), 0);
    
    int ii;
    for( ii = 0; ii < (int)tVec.size(); ii++ )
    {
      double t = tVec[ii];

        //Get the s-u surface at this t value.
      int jj, kk;
      for( kk = 0; kk <= uDegree; kk++ )
      {
        int uOffset = kk*(sDegree+1)*(tDegree+1);
        for( jj = 0; jj <= sDegree; jj++ )
        {
          int tt;
          int suOffset = uOffset+jj;
          for( tt = 0; tt <= tDegree; tt++ )
          {
            int offset = cpDim*(suOffset + tt*(sDegree+1));
            int mm;
            for( mm = 0; mm < cpDim; mm++ )
                locPts.push_back( ctrlPts[offset+mm] );
          }

          if( 0 != deCasteljauCurve( locPts, cpDim, t, pt ) )
          {
            std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 3" << std::endl;
            return 1;
          }
          tempPts.insert( tempPts.end(), pt.begin(), pt.end() );
          
          pt.clear();
          locPts.clear();
        }
      }

        //Now get the points from the surface.
      if( 0 != deCasteljauSurface( tempPts, sDegree, uDegree, cpDim, sVec, uVec, surfPts ) )
      {
        std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 4" << std::endl;
        return 1;
      }

      std::vector<double>::iterator s_iter = surfPts.begin();
      int tOffset = ii*sVec.size();
      for( kk = 0; kk < (int)uVec.size(); kk++ )
      {
        int uOffset = kk*tVec.size()*sVec.size();
        for( jj = 0; jj < (int)sVec.size(); jj++ )
        {
          int offset = cpDim*(uOffset+tOffset+jj);
          int mm;
          for( mm = 0; mm < cpDim; mm++ )
          {
            pts[offset+mm] = *s_iter;
            s_iter++;
          }
        }
      }
      
      tempPts.clear();
      surfPts.clear();
    }  
  }
  else
  {    
    int ii;
    for( ii = 0; ii < (int)uVec.size(); ii++ )
    {
      double u = uVec[ii];

        //Get the s-t surface at this u value.
      int jj, kk;
      for( kk = 0; kk <= tDegree; kk++ )
      {
        int tOffset = kk*(sDegree+1);
        for( jj = 0; jj <= sDegree; jj++ )
        {
          int uu;
          int stOffset = tOffset+jj;
          for( uu = 0; uu <= uDegree; uu++ )
          {
            int offset = cpDim*(stOffset + uu*(sDegree+1)*(tDegree+1));
            int mm;
            for( mm = 0; mm < cpDim; mm++ )
                locPts.push_back( ctrlPts[offset+mm] );
          }

          if( 0 != deCasteljauCurve( locPts, cpDim, u, pt ) )
          {
            std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 5" << std::endl;
            return 1;
          }
          tempPts.insert( tempPts.end(), pt.begin(), pt.end() );
          
          pt.clear();
          locPts.clear();
        }
      }

        //Now get the points from the surface.
      if( 0 != deCasteljauSurface( tempPts, sDegree, tDegree, cpDim, sVec, tVec, surfPts ) )
      {
        std::cerr << "ERROR: in NURBS_T::deCasteljauVolume() -- 5" << std::endl;
        return 1;
      }

      pts.insert( pts.end(), surfPts.begin(), surfPts.end() );
      
      tempPts.clear();
      surfPts.clear();
    }    
  }
  
  return 0;
}

int NURBS_T::deCasteljauSurface(const std::vector<double>& ctrlPts,
                                const int sDegree,
                                const int tDegree,
                                const int Dim,
                                const std::vector<double>& sVec,
                                const std::vector<double>& tVec,
                                std::vector<double>& pts)
{
  pts.clear();
  
  std::vector<double> tempPts;
  std::vector<double> pt;
  
  const int cpDim = Dim;
  if( sDegree <= tDegree )
  {
    pts.resize(cpDim*sVec.size()*tVec.size(), 0);
    std::vector<double>::iterator p_iter;
    
      //Loop over each s value
    int ii;
    for( ii = 0; ii < (int)sVec.size(); ii++ )
    {
      double s = sVec[ii];
      
        //Find column points first.
      std::vector<double> locPts;
      int jj;
      std::vector<double>::const_iterator iter = ctrlPts.begin();
      for( jj = 0; jj <= tDegree; jj++ )
      {
          //Get the jj row of points.
        tempPts.clear();
        tempPts.insert( tempPts.end(), iter, iter+cpDim*(sDegree+1) );
        iter += cpDim*(sDegree+1);

        if( 0 != deCasteljauCurve( tempPts, Dim, s, pt ) )
        {
          std::cerr << "ERROR: in NURBS_T::deCasteljauSurface()" << std::endl;
          return 1;
        }
        locPts.insert( locPts.end(), pt.begin(), pt.end() );

        pt.clear();
      }

        //Now calculate each t point along this curve.
      for( jj = 0; jj < (int)tVec.size(); jj++ )
      {
        double t = tVec[jj];
        if( 0 != deCasteljauCurve( locPts, Dim, t, pt ) )
        {
          std::cerr << "ERROR: in NURBS_T::deCasteljauSurface() -- 1" << std::endl;
          return 1;
        }

        p_iter = pts.begin() + cpDim*(jj*sVec.size()+ii);
        int kk;
        for( kk = 0; kk < cpDim; kk++ )
            *(p_iter + kk ) = pt[kk];
        pt.clear();
      }
    }
  }
  else
  {
      //Loop over each t value first.
    int ii;
    for( ii = 0; ii < (int)tVec.size(); ii++ )
    {
      double t = tVec[ii];
      
        //Find row points first.
      std::vector<double> locPts;
      std::vector<double>::const_iterator iter;
    
      int jj;
      for( jj = 0; jj <= sDegree; jj++ )
      {
        tempPts.clear();

          //Get the jj column of points from the from the kk row. 
        int kk;
        for( kk = 0; kk <= tDegree; kk++ )
        {
          iter = ctrlPts.begin() + cpDim*(kk*(sDegree+1)+jj);
          int mm;
          for( mm = 0; mm < cpDim; mm++ )
          {
            tempPts.push_back( *iter );
            iter++;
          }
        }

        if( 0 != deCasteljauCurve( tempPts, Dim, t, pt ) )
        {
          std::cerr << "ERROR: in NURBS_T::deCasteljauSurface()" << std::endl;
          return 1;
        }
        locPts.insert( locPts.end(), pt.begin(), pt.end() );

        pt.clear();
      }

        //Now calculate each s point long the curve.
      for( jj = 0; jj < (int)sVec.size(); jj++ )
      {
        double s = sVec[jj];
        if( 0 != deCasteljauCurve( locPts, Dim, s, pt ) )
        {
          std::cerr << "ERROR: in NURBS_T::deCasteljauSurface() -- 2" << std::endl;
          return 1;
        }

        pts.insert( pts.end(), pt.begin(), pt.end() );
        pt.clear();
      }
    }
  }
  
  return 0;
}

int NURBS_T::decomposeCurve(const std::vector<double>& knotVec,
                            const int cDegree,
                            const std::vector<double>& ctrlPts,
                            const int Dim,
                            std::vector<double>& bCtrlPts,
                            std::vector<BezierElem*>& bSegments,
                            const bool calcbCtrlPts)
{
  bCtrlPts.clear();
  bSegments.clear();
  
  const int cpDim = Dim;
    //Short circuit if curve is made up of only one segement
    //to begin with.
  if( (int)ctrlPts.size()/cpDim == cDegree + 1)
  {

    BezierElem* elem = new BezierElem( cDegree );
    bSegments.push_back( elem );
    
    int ii;
    for( ii = 0; ii < (int)ctrlPts.size()/cpDim; ii++ )
    {
      elem->gPtsMap.push_back(ii);
      elem->lPtsMap.push_back(ii);
    }

    if( calcbCtrlPts )
    {
      bCtrlPts.insert( bCtrlPts.end(), ctrlPts.begin(), ctrlPts.end() );
    }
    
    return 0;
  }
  
  int a = cDegree;
  int b = a + 1;
  int numSegs = 0;

  int ii, mm;
  
    //Initialize first segment.
  BezierElem* curr_elem = new BezierElem( cDegree );
  BezierElem* next_elem = NULL;
  bSegments.push_back(curr_elem);
  for( ii = 0; ii <= cDegree; ii++ )
  {
    if( calcbCtrlPts )
    {
      for( mm = 0; mm < cpDim; mm++ )
          bCtrlPts.push_back( ctrlPts[cpDim*ii+mm] );
    }
    
    curr_elem->gPtsMap.push_back(ii);
    curr_elem->lPtsMap.push_back(ii);
  }
  
  while( b < (int)knotVec.size()-1 )
  {
    const int offset = numSegs*(cDegree);
    ii = b;

      //Count knot multiplicities.
    while( b < (int)knotVec.size()-1 && knotVec[b+1] == knotVec[b] )
        b++;

    if( b < (int)knotVec.size()-1 )
    {
        //Allocate space for the next segment.
      if( calcbCtrlPts )
          bCtrlPts.insert(bCtrlPts.end(), cDegree*cpDim, 0.0 );

      next_elem = new BezierElem(cDegree);
      bSegments.push_back(next_elem);
    }
    
    const int mult = b-ii+1;
    if( mult < cDegree )
    {
        //Insert new knots to increase multiplicity to cDegree.
      double numer = knotVec[b] - knotVec[a];

        //Compute and store alphas.
      std::vector<double> alphas;
      alphas.resize( cDegree - 1, 0.0 );
      int jj;
      for( jj = cDegree; jj > mult; jj-- )
          alphas[jj-mult-1] = numer/ (knotVec[a+jj]-knotVec[a] );

        //Insert knot cDegree-mult times.
      const int r = cDegree-mult;
      for( jj = 1; jj <= r; jj++ )
      {
        int save = r - jj;
        const int s = mult+jj;  //Number of new points.
        int kk;
        for( kk = cDegree; kk >= s; kk-- )
        {
          double alpha = alphas[kk-s];

          if( calcbCtrlPts )
          {
            for( mm = 0; mm < cpDim; mm++ )
                bCtrlPts[cpDim*(offset+kk)+mm] = alpha*bCtrlPts[cpDim*(offset+kk)+mm]+
                    (1.0-alpha)*bCtrlPts[cpDim*(offset+kk-1)+mm];
          }
          
            //Update matrix coefficients.
          for( mm = 0; mm <= cDegree; mm++ )
          {
            curr_elem->coefficient( mm, kk,
                                    alpha*curr_elem->coefficient(mm, kk)
                                    + (1.0-alpha)*curr_elem->coefficient(mm, kk-1) );
          }
          
        }

          //Update control points of next segment.
        if( b < (int)knotVec.size()-1 )
        {
          if( save > 0 && calcbCtrlPts )
          {
            for( mm = 0; mm < cpDim; mm++ )
                bCtrlPts[cpDim*(offset+cDegree+save)+mm] =
                    bCtrlPts[cpDim*(offset+cDegree)+mm];
          }

            // Update coefficient Matrix
          for( mm = 0; mm <= jj; mm++ )
              next_elem->coefficient( mm+save, save,
                                      curr_elem->coefficient( cDegree-jj+mm, cDegree ) );
        }
      }
    }

      //Segment complete
    numSegs++;

    if( b < (int)knotVec.size()-1 )
    {
        //Initialize for next segment.
      int lastgI = curr_elem->gPtsMap.back() - cDegree + mult;
      int lastlI = curr_elem->lPtsMap.back();
      
      if( calcbCtrlPts )
      {
        for( ii = cDegree-mult; ii <= cDegree; ii++ )
        {
          for( mm = 0; mm < cpDim; mm++ )
              bCtrlPts[cpDim*(offset+cDegree+ii)+mm] = ctrlPts[cpDim*(b-cDegree+ii)+mm];
        }
      }
      
      for( ii = 0; ii <= cDegree; ii++ )
      {
        next_elem->gPtsMap.push_back(lastgI+ii);
        next_elem->lPtsMap.push_back(lastlI+ii);
      }
      
      a = b;
      b++;

      curr_elem = next_elem;
    }
  }
  
  return 0;
}

int NURBS_T::decomposeSurface(const std::vector<double>& sKnots,
                              const int sDegree,
                              const std::vector<double>& tKnots,
                              const int tDegree,
                              const std::vector<double>& ctrlPts,
                              const int Dim,
                              std::vector<double>& bCtrlPts,
                              std::vector<BezierElem*>& bSegments,
                              const bool buildFullOps)
{
  bCtrlPts.clear();
  bSegments.clear();
  
    //Do a decomposition of each knot vector individually.  This will give
    //the 1D coefficient matrices that can then be used to get the coefficient
    //matrices for the surface Beziers.  Once these coefficients are obtained
    //they can be used to calculated the control points.  Not sure if this is
    //the most efficient calculation but it is the most efficient use of
    //existing code.
  std::vector<BezierElem*> sSegments;
  std::vector<double> sCtrlPts;
  std::vector<double> sbCtrlPts;
  std::vector<BezierElem*> tSegments;
  std::vector<double> tCtrlPts;
  std::vector<double> tbCtrlPts;

  const int numSPts = sKnots.size() - sDegree - 1;
  const int numTPts = tKnots.size() - tDegree - 1;
  const int cpDim = Dim;

    //Get the first row of points in the s direction and
    //first column in the t direction.
  int ii, jj;
  for( ii = 0; ii < numSPts; ii++ )
      for( jj = 0; jj < cpDim; jj++ )
          sCtrlPts.push_back( ctrlPts[cpDim*ii+jj] );
  for( ii = 0; ii < numTPts; ii++ )
      for( jj = 0; jj < cpDim; jj++ )
          tCtrlPts.push_back( ctrlPts[cpDim*(ii*numSPts)+jj] );
  
  int rc = 0; //Return condition
  
  do //Loop used to control memory de-allocation.
  {
    if( 0 != decomposeCurve( sKnots, sDegree, sCtrlPts, Dim,
                             sbCtrlPts, sSegments ) )
    {
      std::cerr << "Unable to decompose surface in s"
                << "direction in NURBS_T::decomposeSurfac()" << std::endl;
      rc = 1;
      break;
    }
    
    if( 0 != decomposeCurve( tKnots, tDegree, tCtrlPts, Dim,
                             tbCtrlPts, tSegments ) )
    {
      std::cerr << "Unable to decompose surface in t"
                << "direction in NURBS_T::decomposeSurfac()" << std::endl;
      rc = 1;
      break;
    }

      //Add the first row of control points to the final list.
    bCtrlPts.insert( bCtrlPts.end(), sbCtrlPts.begin(), sbCtrlPts.end() );

    const int numSBPts = sSegments.size()*sDegree+1;
    const int numTBPts = tSegments.size()*tDegree+1;
    
      //Resize the bCtrlPts vector so that array indexing will point to
      //valid slots.
    bCtrlPts.resize( cpDim*( numSBPts*numTBPts ) , 0.0 );

      //Now add the first column of control points to the final list.
    for( ii = 1; ii < numTBPts; ii++ )
        for( jj = 0; jj < cpDim; jj++ )
            bCtrlPts[cpDim*(ii*numSBPts)+jj] = tbCtrlPts[cpDim*ii+jj];
  
    int sIndex, tIndex;
    for( tIndex = 0; tIndex < (int)tSegments.size(); tIndex++ )
    {
      for( sIndex = 0; sIndex < (int)sSegments.size(); sIndex++ )
      {
          //Create the new segment.
        BezierElem* elem = new BezierElem( sDegree, tDegree, 0, buildFullOps );
        bSegments.push_back( elem );

        if( false == buildFullOps )
        {
          int row, col;
          for( row = 0; row <= sDegree; row++ )
            for( col = 0; col <= sDegree; col++ )
              elem->set_S_MatrixCoeff( row, col, sSegments[sIndex]->coefficient( row, col ) );
          
          for( row = 0; row <= tDegree; row++ )
            for( col = 0; col <= tDegree; col++ )
              elem->set_T_MatrixCoeff( row, col, tSegments[tIndex]->coefficient( row, col ) );
        }
        else
        {
          //Tensor Product to create new coefficient matrix (T * S).
          int tRow, tCol;
          for( tRow = 0; tRow <= tDegree; tRow++ )
          {
            for( tCol = 0; tCol <= tDegree; tCol++ )
            {
              double t = tSegments[tIndex]->coefficient( tRow, tCol );
              int rowOffset = tRow * (sDegree+1);
              int colOffset = tCol * (sDegree+1);
              
              int sRow, sCol;
              for( sRow = 0; sRow <= sDegree; sRow++ )
              {
                for( sCol = 0; sCol <= sDegree; sCol++ )
                {
                  elem->coefficient( sRow + rowOffset, sCol + colOffset,
                                    t * sSegments[sIndex]->coefficient( sRow, sCol ) );
                }
              }
            }
          }
        }  
          //Create the map to local and global control points.
        int ss, tt;
        for( tt = 0; tt <= tDegree; tt++ )
        {
          for( ss = 0; ss <= sDegree; ss++ )
          {
            elem->gPtsMap.push_back( tSegments[tIndex]->gPtsMap[tt]*numSPts
                                     + sSegments[sIndex]->gPtsMap[ss] );

            elem->lPtsMap.push_back( (tIndex*tDegree+tt)*numSBPts + sIndex*sDegree+ss );
          }
        }                             
          
          //Use map and coefficients to create new control points.  Note: min
          //s and t points have already been create by last pass.
        int cRow, cCol;
        int rowOffset = 1;
        int colOffset = 1;
        for( cRow = sDegree + 2; cRow < (sDegree+1)*(tDegree+1); cRow++ )
        {
          int bPtOffset = ( (tIndex*tDegree+rowOffset)*numSBPts
                            + sIndex*sDegree + colOffset );
            
          int  mm;
          for( cCol = 0; cCol < (sDegree+1)*(tDegree+1); cCol++ )
          {
              //Note: multiplying by matrix transpose.
            double coeff = elem->coefficient(cCol, cRow);
            if( 0.0 != coeff )
            {
              for( mm = 0; mm < cpDim; mm++ )
                  bCtrlPts[cpDim*(bPtOffset)+mm] += coeff * ctrlPts[cpDim*elem->gPtsMap[cCol]+mm];
            }
          }

            //Skip row containing point at min t.
          if( sDegree == cRow % (sDegree+1) )
          {
            cRow++;
            rowOffset++;
            colOffset = 1;
          }
          else
              colOffset++;
        }
      }
    }
    
  }while( 0 );  //Only one pass.
  
    //Delete memory allocated by curve decomposition.
  while( !sSegments.empty() )
  {
    delete sSegments.back();
    sSegments.pop_back();
  }
  while( !tSegments.empty() )
  {
    delete tSegments.back();
    tSegments.pop_back();
  }
      
  return rc;
}

int NURBS_T::decomposeVolume(const std::vector<double>& sKnots,
                             const int sDegree,
                             const std::vector<double>& tKnots,
                             const int tDegree,
                             const std::vector<double>& uKnots,
                             const int uDegree,
                             const std::vector<double>& ctrlPts,
                             const int Dim,
                             std::vector<double>& bCtrlPts,
                             std::vector<BezierElem*>& bSegments,
                             const bool buildFullOps)
{
  bCtrlPts.clear();
  bSegments.clear();

    //Do a decomposition of each knot vector individually.  This will give
    //the 1D coefficient matrices that can then be used to get the coefficient
    //matrices for the surface Beziers.  Once these coefficients are obtained
    //they can be used to calculated the control points.  Not sure if this is
    //the most efficient calculation but it is the most efficient use of
    //existing code.
  std::vector<BezierElem*> sSegments;
  std::vector<double> sCtrlPts;
  std::vector<double> sbCtrlPts;
  std::vector<BezierElem*> tSegments;
  std::vector<double> tCtrlPts;
  std::vector<double> tbCtrlPts;
  std::vector<BezierElem*> uSegments;
  std::vector<double> uCtrlPts;
  std::vector<double> ubCtrlPts;

  const int numSPts = sKnots.size() - sDegree - 1;
  const int numTPts = tKnots.size() - tDegree - 1;
  const int numUPts = uKnots.size() - uDegree - 1;
  const int cpDim = Dim;

    //Get the initial points in the s, t, and u direction.
  int ii, jj;
  std::vector<double>::const_iterator c_iter;
  for( ii = 0; ii < numSPts; ii++ )
  {
    c_iter = ctrlPts.begin() + cpDim*ii;
    for( jj = 0; jj < cpDim; jj++ )
    {
      sCtrlPts.push_back( *c_iter );
      c_iter++;
    }
  }
  for( ii = 0; ii < numTPts; ii++ )
  {
    c_iter = ctrlPts.begin() + cpDim*(ii*numSPts);
    for( jj = 0; jj < cpDim; jj++ )
    {
      tCtrlPts.push_back( *c_iter );
      c_iter++;
    }
  }
  for( ii = 0; ii < numUPts; ii++ )
  {
    c_iter = ctrlPts.begin() + cpDim*(ii*numSPts*numTPts);
    for( jj = 0; jj < cpDim; jj++ )
    {
      uCtrlPts.push_back( *c_iter );
      c_iter++;
    }
  }
  
  int rc = 0; //Return condition
  
  do //Loop used to control memory de-allocation.
  {
    if( 0 != decomposeCurve( sKnots, sDegree, sCtrlPts, cpDim,
                             sbCtrlPts, sSegments ) )
    {
      std::cerr << "Unable to decompose surface in s"
                << "direction in NURBS_T::decomposeSurfac()" << std::endl;
      rc = 1;
      break;
    }
    
    if( 0 != decomposeCurve( tKnots, tDegree, tCtrlPts, cpDim,
                             tbCtrlPts, tSegments ) )
    {
      std::cerr << "Unable to decompose surface in t"
                << "direction in NURBS_T::decomposeSurfac()" << std::endl;
      rc = 1;
      break;
    }

    if( 0 != decomposeCurve( uKnots, uDegree, uCtrlPts, cpDim,
                             ubCtrlPts, uSegments ) )
    {
      std::cerr << "Unable to decompose surface in u"
                << "direction in NURBS_T::decomposeSurfac()" << std::endl;
      rc = 1;
      break;
    }

      //Add the first row of control points to the final list.
    bCtrlPts.insert( bCtrlPts.end(), sbCtrlPts.begin(), sbCtrlPts.end() );

    const int numSBPts = sSegments.size()*sDegree+1;
    const int numTBPts = tSegments.size()*tDegree+1;
    const int numUBPts = uSegments.size()*uDegree+1;
    
      //Resize the bCtrlPts vector so that array indexing will point to
      //valid slots.
    bCtrlPts.resize( cpDim*( numSBPts*numTBPts*numUBPts ) , 0.0 );

    std::vector<double>::iterator r_iter;  //RHS
    std::vector<double>::iterator l_iter;  //LHS

      //Now add the first column of control points to the final list.    
    for( ii = 1; ii < numTBPts; ii++ )
    {
      l_iter = bCtrlPts.begin() + cpDim*(ii*numSBPts);
      r_iter = tbCtrlPts.begin() + cpDim*ii;
      for( jj = 0; jj < cpDim; jj++ )
          *(l_iter+jj) = *(r_iter+jj);
    }
    
      //Now the first points in u direction.
    for( ii = 1; ii < numUBPts; ii++ )
    {
      l_iter = bCtrlPts.begin() + cpDim*(ii*numSBPts*numTBPts);
      r_iter = ubCtrlPts.begin() + cpDim*ii;
      for( jj = 0; jj < cpDim; jj++ )
            *(l_iter+jj) = *(r_iter+jj);
    }
    
    int sIndex, tIndex, uIndex;
    for( uIndex = 0; uIndex < (int)uSegments.size(); uIndex++ )
    {
      for( tIndex = 0; tIndex < (int)tSegments.size(); tIndex++ )
      {
        for( sIndex = 0; sIndex < (int)sSegments.size(); sIndex++ )
        {
            //Create the new segment.
          BezierElem* elem = new BezierElem( sDegree, tDegree, uDegree, buildFullOps );
          bSegments.push_back( elem );
          
          if( false == buildFullOps )
          {
            int row, col;
            for( row = 0; row <= sDegree; row++ )
              for( col = 0; col <= sDegree; col++ )
                elem->set_S_MatrixCoeff( row, col, sSegments[sIndex]->coefficient( row, col ) );
            
            for( row = 0; row <= tDegree; row++ )
              for( col = 0; col <= tDegree; col++ )
                elem->set_T_MatrixCoeff( row, col, tSegments[tIndex]->coefficient( row, col ) );
            
            for( row = 0; row <= uDegree; row++ )
              for( col = 0; col <= uDegree; col++ )
                elem->set_U_MatrixCoeff( row, col, uSegments[uIndex]->coefficient( row, col ) );
          }
          else
          {
            //Tensor Product to create new coefficient matrix: U * (T * S).
            int uRow, uCol;
            for( uRow = 0; uRow <= uDegree; uRow++ )
            {
              for( uCol = 0; uCol <= uDegree; uCol++ )
              {
                double u = uSegments[uIndex]->coefficient( uRow, uCol );
                
                if( 0.0 == u )
                  continue;
                
                int tRow, tCol;
                for( tRow = 0; tRow <= tDegree; tRow++ )
                {
                  for( tCol = 0; tCol <= tDegree; tCol++ )
                  {
                    double t = tSegments[tIndex]->coefficient( tRow, tCol );
                    
                    if( 0.0 == t )
                      continue;
                    
                    int rowOffset = uRow*(sDegree+1)*(tDegree+1) + tRow * (sDegree+1);
                    int colOffset = uCol*(sDegree+1)*(tDegree+1) + tCol * (sDegree+1);
                    
                    int sRow, sCol;
                    for( sRow = 0; sRow <= sDegree; sRow++ )
                    {
                      for( sCol = 0; sCol <= sDegree; sCol++ )
                      {
                        elem->coefficient( sRow + rowOffset, sCol + colOffset,
                                          u * t * sSegments[sIndex]->coefficient( sRow, sCol ) );
                      }
                    }
                  }
                }
              }
            }
          }
          //Create the map to local and global control points.
          int ss, tt, uu;
          for( uu = 0; uu <= uDegree; uu++ )
          {
            for( tt = 0; tt <= tDegree; tt++ )
            {
              for( ss = 0; ss <= sDegree; ss++ )
              {
                elem->gPtsMap.push_back( uSegments[uIndex]->gPtsMap[uu]*numSPts*numTPts
                                         + tSegments[tIndex]->gPtsMap[tt]*numSPts
                                         + sSegments[sIndex]->gPtsMap[ss] );

                elem->lPtsMap.push_back( (uIndex*uDegree+uu)*numSBPts*numTBPts
                                         + (tIndex*tDegree+tt)*numSBPts
                                         + sIndex*sDegree+ss );
              }
            }                             
          }
          
            //Use map and coefficients to create new control points.
          const int matSize = (sDegree+1)*(tDegree+1)*(uDegree+1);
          
            //u = 0 plane
          int cRow, cCol;
          int rowOffset = 0, colOffset = 0;
          if( 0 == uIndex )
          {
            rowOffset = 1;
            colOffset = 1;
            for( cRow = sDegree + 2; cRow < (sDegree+1)*(tDegree+1); cRow++ )
            {
              int bPtOffset = ( (tIndex*tDegree+rowOffset)*numSBPts
                                + sIndex*sDegree + colOffset );
            
              int  mm;
              for( cCol = 0; cCol < matSize; cCol++ )
              {
                  //Note: multiplying by matrix transpose.
                double coeff = elem->coefficient(cCol, cRow);
                if( 0.0 != coeff )
                {
                  l_iter = bCtrlPts.begin() + cpDim*bPtOffset;
                  c_iter = ctrlPts.begin() + cpDim*elem->gPtsMap[cCol];
                  for( mm = 0; mm < cpDim; mm++ )
                      *(l_iter+mm) += coeff * (*(c_iter+mm));
                }
              }

                //Skip row containing point at min t.
              if( sDegree == cRow % (sDegree+1) )
              {
                cRow++;
                rowOffset++;
                colOffset = 1;
              }
              else
                  colOffset++;
            }
          }
          else
              cRow = (sDegree+1)*(tDegree+1)+1;

          cRow--;
          
            //u != 0 planes, cRow starts at the first point of the second layer.
            //(first point is on the s=t=0 curve).
          int layerOffset = 0;
          for( ; cRow < matSize; cRow++ )
          {
            if( 0 == cRow % ((sDegree+1)*(tDegree+1)) )
            {
                //This code should be hit on the first pass through of 
                //each layer.
              cRow++;
              
              if( 0 != tIndex )
              {
                  //Skip the "bottom" points.
                cRow += sDegree;
                rowOffset = 1;
                colOffset = 0;
              }
              else
              {
                rowOffset = 0;
                colOffset = 1;
              }

              layerOffset++;
            }
            
            if( 0 != sIndex && 0 == cRow % (sDegree+1) )
            {
                //Skip row for point on s-u plane.
              cRow++;
              colOffset++;
            }
            
            int bPtOffset = ( (uIndex*uDegree+layerOffset)*numTBPts*numSBPts
                              + (tIndex*tDegree+rowOffset)*numSBPts
                              + sIndex*sDegree + colOffset );
            
            int mm;
            for( cCol = 0; cCol < matSize; cCol++ )
            {
                //Note: multiplying by matrix transpose.
              double coeff = elem->coefficient(cCol, cRow);
              if( 0.0 != coeff )
              {
                l_iter = bCtrlPts.begin() + cpDim*bPtOffset;
                c_iter = ctrlPts.begin() + cpDim*elem->gPtsMap[cCol];
                for( mm = 0; mm < cpDim; mm++ )
                    *(l_iter+mm) += coeff * (*(c_iter+mm));
              }
            }

            if( sDegree == cRow % (sDegree+1) )
            {
              rowOffset++;
              colOffset = 0;
            }
            else
                colOffset++;
          }
        }
      }
    }
  }while( 0 );  //Only one pass.

  
   //Delete memory allocated by curve decomposition.
  while( !sSegments.empty() )
  {
    delete sSegments.back();
    sSegments.pop_back();
  }
  while( !tSegments.empty() )
  {
    delete tSegments.back();
    tSegments.pop_back();
  }
  while( !uSegments.empty() )
  {
    delete uSegments.back();
    uSegments.pop_back();
  }
  
  return rc;
}

int NURBS_T::transformPoints(const std::vector<double>& inPts,
                             const std::vector<BezierElem*>& bSegments,
                             const int Dim,
                             std::vector<double>& outPts)
{
  outPts.clear();
  
    //Now transform the points to the Bezier segments using
    //the coefficient matrices.
  std::map<const int, double*> pointMap;
  int ii;
  for( ii = 0; ii < (int)bSegments.size(); ii++ )
  {
    NURBS_T::BezierElem* elem = bSegments[ii];
    const int matSize = (elem->sDegree()+1)*(elem->tDegree()+1)*(elem->uDegree()+1);

      //Loop through the transpose of the coefficient matrix.
    int row, col;
    for( row = 0; row < matSize; row++ )
    {
      int bPoint = elem->lPtsMap[row];
      if( pointMap.find(bPoint) == pointMap.end() )
      {
          //Create the new point.
        double* pt = new double[4];
        pointMap[bPoint] = pt;
        pt[0] = pt[1] = pt[2] = pt[3] = 0;
        for( col = 0; col < matSize; col++ )
        {
          double coeff = elem->coefficient(col, row);
          if( coeff != 0 )
          {
            std::vector<double>::const_iterator d_iter = inPts.begin() + 4*elem->gPtsMap[col];

            int mm;
            for( mm = 0; mm < 4; mm++ )
                pt[mm] += coeff * *(d_iter+mm);
          }
        }
      }
    }
  }

    //Now copy to point into the return vector.
  for( ii = 0; ii < (int)pointMap.size(); ii++ )
  {
    assert( pointMap.find( ii ) != pointMap.end() );
    
    double* pt = pointMap[ii];
    int mm;
    for( mm = 0; mm < 4; mm++ )
        outPts.push_back( pt[mm] );

    delete [] pt;
  }
  
  return 0;
}

int NURBS_T::basisFunctions(const std::vector<double>& knotVec,
                            const int cDegree,
                            const double u,
                            std::vector<double>& N)
{
  N.clear();
  N.resize( cDegree+1, 0.0);
  
  int span = findSpan(knotVec, cDegree, u);
  if( span <  0 )
    return span;
  
  double* left = new double[cDegree+1];
  double* right = new double[cDegree+1];
  
  int jj;
  N[0] = 1.0;
  for( jj = 1; jj <= cDegree; jj++ )
  {
    left[jj] = u - knotVec[span+1-jj];
    right[jj] = knotVec[span+jj] - u;
    
    double saved = 0.0;
    int r;
    for( r = 0; r < jj; r++ )
    {
      double temp = N[r]/(right[r+1]+left[jj-r]);
      N[r] = saved + right[r+1]*temp;
      saved = left[jj-r]*temp;
    }
    N[jj] = saved;
  }
  
  delete [] left;
  delete [] right;
  return span;
}

int NURBS_T::basisFuncsAndDers(const std::vector<double>& knotVec,
                               const int cDegree,
                               const double u,
                               const int numDers,
                               std::vector<double>& ders)
{
  ders.clear();
  
  //Find the span in which u lies.  Return if
  //u is not within the limits of the knot vector.
  int span = findSpan( knotVec, cDegree, u);
  if( span < 0 )
    return span;
  
  //Don't calculate zero derivatives.
  int nDers = numDers;
  if( nDers > cDegree )
    nDers = cDegree;
  
  ders.resize( (cDegree+1)*(nDers+1), 0.0 );
  
  int j, r;
  
  double** ndu = new double*[cDegree+1];
  for( j = 0; j <= cDegree; j++ )
    ndu[j] = new double[cDegree+1];
  double* left = new double[cDegree+1];
  double* right = new double[cDegree+1];
  
  //Calculate basis functions and knot differences.
  ndu[0][0] = 1.0;
  for( j = 1; j <= cDegree; j++ )
  {
    left[j] = u - knotVec[span+1-j];
    right[j] = knotVec[span+j] - u;
    
    double saved = 0.0;
    for( r = 0; r < j; r++ )
    {
      ndu[j][r] = right[r+1]+left[j-r];
      double temp = ndu[r][j-1]/ndu[j][r];
      
      ndu[r][j] = saved + right[r+1]*temp;
      saved = left[j-r]*temp;
    }
    ndu[j][j] = saved;
  }
  //Load the basis functions.
  for( j = 0; j <= cDegree; j++ )
    ders[j*(nDers+1)] = ndu[j][cDegree];
  
  double* a[2];
  a[0] = new double[cDegree+1];
  a[1] = new double[cDegree+1];
  
  //Compute the derivatives.
  int k;
  for( r = 0; r <= cDegree; r++ )
  {
    int s1 = 0;
    int s2 = 1;
    
    a[0][0] = 1.0;
    
    //Compute the kth derivatives.
    for( k = 1; k <= nDers; k++ )
    {
      double d = 0.0;
      int rk = r-k;
      int pk = cDegree-k;
      if( r >= k )
      {
        a[s2][0] = a[s1][0]/ndu[pk+1][rk];
        d = a[s2][0]*ndu[rk][pk];
      }
      
      int j1, j2;
      if( rk >= -1 )
        j1 = 1;
      else
        j1 = -rk;
      
      if( r-1 <= pk )
        j2 = k-1;
      else
        j2 = cDegree - r;
      
      for( j = j1; j <= j2; j++ )
      {
        a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pk+1][rk+j];
        d += a[s2][j]*ndu[rk+j][pk];
      }
      
      if( r <= pk )
      {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
        d += a[s2][k]*ndu[r][pk];
      }
      
      ders[r*(nDers+1)+k] = d;
      
      //switch rows.
      j = s1; 
      s1 = s2;
      s2 = j;
    }
  }
  
  //Multiply through be the correct factors.
  r = cDegree;
  for( k = 1; k <= nDers; k++ )
  {
    for( j = 0; j <= cDegree; j++ )
      ders[j*(nDers+1)+k] *= r;
    
    r*= (cDegree-k);
  }
  
  
  //Clean up memory
  for( j = 0; j <= cDegree; j++ )
    delete [] ndu[j];
  delete [] ndu;
  delete [] a[0];
  delete [] a[1];
  delete [] left;
  delete [] right;
  
  return span;
}

int NURBS_T::computeGrevilleAbscissae(const std::vector<double>& knotVec,
                                      const int cDegree,
                                      std::vector<double>& pts)
{
  //Will compute knotVec.size() - cDegree - 1 values.
  int stop = knotVec.size() - cDegree;
  int ii;
  for( ii = 1; ii < stop; ii++ )
  {
    double sum = 0.0;
    int jj;
    for( jj = 0; jj < cDegree; jj++ )
      sum += knotVec[ii+jj];
    sum = sum / double( cDegree );
    
    pts.push_back( sum );
  }
  assert( pts.size() == (knotVec.size() - cDegree-1 ) );
  
  return 0;
}

int NURBS_T::degreeElevateCurve(std::vector<double>& oldKnotVec,
                                const int cDegree,
                                const int Dim,
                                std::vector<double>& oldCtrlPts,
                                const int addDegrees)
{
  std::vector<double> newKnotVec;
  std::vector<double> newPts;
  
  if( 0 != degreeElevateCurve( oldKnotVec, cDegree, Dim,
                               oldCtrlPts, addDegrees,
                               newKnotVec, newPts ) )
    return 1;
  
  oldKnotVec.clear();
  oldCtrlPts.clear();
  oldKnotVec.insert( oldKnotVec.end(), newKnotVec.begin(), newKnotVec.end() );
  oldCtrlPts.insert( oldCtrlPts.end(), newPts.begin(), newPts.end() );
  
  return 0;
}

int NURBS_T::degreeElevateCurve(const std::vector<double>& oldKnotVec,
                                const int cDegree,
                                const int Dim,
                                const std::vector<double>& oldCtrlPts,
                                const int addDegrees,
                                std::vector<double>& newKnotVec,
                                std::vector<double>& newCtrlPts)
{
  int ii, jj, dd;
  newKnotVec.clear();
  newCtrlPts.clear();
  
  const int m = oldKnotVec.size() - 1;
  const int ph = cDegree + addDegrees;
  const int ph2 = ph/2;
  
  //Compute Bezier degree elevation coefficients.
  double** bezalfs = new double*[ph+1];
  for( ii = 0; ii <= ph; ii++ )
    bezalfs[ii] = new double[cDegree+1];
  
  bezalfs[0][0] = bezalfs[ph][cDegree] = 1.0;
  for( ii = 1; ii <= ph2; ii++ )
  {
    //Note: If this is too expensive then cache binomial coefficients 
    //so that dublicates are not recomputed.
    double inv = 1.0 / double( binomialCoefficient( ph, ii ) );
    int mpi = cDegree < ii ? cDegree : ii;
    
    jj = 0 > ii-addDegrees ? 0 : ii-addDegrees;
    for( ; jj <= mpi; jj++ )
      bezalfs[ii][jj] = inv*binomialCoefficient(cDegree, jj ) 
        * binomialCoefficient( addDegrees, ii - jj );
  }
    
  for( ii = ph2+1; ii <= ph-1; ii++ )
  {
    int mpi = cDegree < ii ? cDegree : ii;
    jj = 0 > ii - addDegrees ? 0 : ii - addDegrees;
    for( ; jj <= mpi; jj++ )
      bezalfs[ii][jj] = bezalfs[ph-ii][cDegree-jj];
  }
  
  int mh = ph;
  int kind = ph+1;
  int r = -1;
  int a = cDegree;
  int b = cDegree+1;
  int cind = 1;
  
  //Initialize first points.
  double ua = oldKnotVec[0];
  for( ii = 0; ii < Dim; ii++ )
    newCtrlPts.push_back( oldCtrlPts[ii] );
  for( ii = 0; ii <= ph; ii++ )
    newKnotVec.push_back( ua );
  
  //Initialize first segment.
  double* bpts = new double[Dim*(cDegree+1)];
  for( ii = 0; ii <= cDegree; ii++ )
    for( jj = 0; jj < Dim; jj++ )
      bpts[ii*Dim+jj] = oldCtrlPts[ii*Dim+jj];
  
  //Big loop thru knot vector
  double* alphas = new double[cDegree-1];
  double* ebpts = new double[(cDegree+addDegrees+1)*Dim];
  double* Nextbpts = new double[(cDegree-1)*Dim];
  
  while( b < m )
  {
    ii = b;
    while( b < m && oldKnotVec[b] == oldKnotVec[b+1])
      b++;
    
    int mul = b-ii+1;
    mh = mh+mul+addDegrees;
    double ub = oldKnotVec[b];
    
    int oldr = r;
    r = cDegree - mul;
    
    //Insert knot u(b) r times.
    int lbz, rbz;
    if( oldr > 0 )
      lbz = (oldr+2)/2;
    else lbz = 1;
    
    if( r > 0 )
      rbz = ph - (r+1)/2;
    else
      rbz = ph;
    
    if( r > 0 )
    {
      //Insert knot to get Bezier segment.
      double numer = ub - ua;
      int kk;
      for( kk = cDegree; kk > mul; kk-- )
        alphas[kk-mul-1] = numer / (oldKnotVec[a+kk] - ua );
      
      for( jj = 1; jj <= r; jj++ )
      {
        int save = r-jj;
        int s = mul+jj;
        for( kk = cDegree; kk >= s; kk-- )
          for( dd = 0; dd < Dim; dd++ )
            bpts[kk*Dim+dd] = alphas[kk-s]*bpts[kk*Dim+dd] + (1.0-alphas[kk-s])*bpts[(kk-1)*Dim+dd];
        
        for( dd = 0; dd < Dim; dd++ )
          Nextbpts[save*Dim+dd] = bpts[cDegree*Dim+dd];
      }
    }
    
    //Degree elevate Bezier.
    for( ii = lbz; ii <= ph; ii++ )
    {
      for( dd = 0; dd < Dim; dd++ )
        ebpts[ii*Dim+dd] = 0.0;
      int mpi = cDegree < ii ? cDegree : ii;
      jj = 0 > ii - addDegrees ? 0 : ii - addDegrees;
      for( ; jj <= mpi; jj++ )
        for( dd = 0; dd < Dim; dd++ )
          ebpts[ii*Dim+dd] += bezalfs[ii][jj]*bpts[jj*Dim+dd];
    }
    
    if( oldr > 1 )
    {
      //Must remove knot u=ua oldr times.
      int first = kind - 2;
      int last = kind;
      
      double den = ub - ua;
      double bet = (ub - newKnotVec[kind-1])/den;
      int tr;
      //Knot removal loop
      for( tr = 1; tr < oldr; tr++ )
      {
        ii = first;
        jj = last;
        int kj = jj - kind + 1;
        
        while( jj - ii > tr )
        {
          // Loop and compute the new control points
          // for one removal step
          if( ii < cind )
          {
            double alf = (ub - newKnotVec[ii])/(ua - newKnotVec[ii]);
            for( dd = 0; dd < Dim; dd++ )
              newCtrlPts[ii*Dim+dd] = alf*newCtrlPts[ii*Dim+dd] + (1.0 - alf)*newCtrlPts[(ii-1)*Dim+dd];
          }
          
          if( jj >= lbz )
          {
            if( jj - tr <= kind - ph + oldr )
            {
              double gam = (ub - newKnotVec[jj-tr])/den;
              for( dd = 0; dd < Dim; dd++ )
                ebpts[kj*Dim+dd] = gam*ebpts[kj*Dim+dd]+(1.0-gam)*ebpts[(kj+1)*Dim+dd];
            }
            else
            {
              for( dd = 0; dd < Dim; dd++ )
                ebpts[kj*Dim+dd] = bet*ebpts[kj*Dim+dd] + (1.0-bet)*ebpts[(kj+1)*Dim+dd];
            }
          }
          
          ii++;
          jj--;
          kj--;
        }
        first--;
        last++;
      }
    }
    
    //Load knot ua.
    if( a != cDegree )
    {
      for( ii = 0; ii < ph-oldr; ii++ )
      {
        newKnotVec.push_back( ua );
        kind++;
      }
    }
    
    
    //Load ctrl points into newCtrlPts.
    for( jj = lbz; jj <= rbz; jj++ )
    {
      for( dd = 0; dd < Dim; dd++ )
        newCtrlPts.push_back( ebpts[jj*Dim+dd] );
      cind++;
    }
    
    if( b < m )
    {
      //Setup for next pass.
      for( jj = 0; jj < r; jj++ )
        for( dd = 0; dd < Dim; dd++ )
          bpts[jj*Dim+dd] = Nextbpts[jj*Dim+dd];
      
      for( jj = r; jj <= cDegree; jj++ )
        for( dd = 0; dd < Dim; dd++ )
          bpts[jj*Dim+dd] = oldCtrlPts[(b-cDegree+jj)*Dim+dd];
      
      a = b;
      b++;
      ua = ub;
    }
    else
    {
      //end knot
      for( ii = 0; ii <= ph; ii++ )
        newKnotVec.push_back( ub );
    }
  }
    
  //Cleanup memory
  for( ii = 0; ii <= ph; ii++ )
    delete [] bezalfs[ii];
  delete [] bezalfs;
  
  delete [] bpts;
  delete [] alphas;
  delete [] Nextbpts;
  delete [] ebpts;
  
  return 0;
}

int NURBS_T::binomialCoefficient(const int n, const int m)
{
  if( m > n )
    return 0;
  
  int locM = m;
  //Take advantage of symmetry.
  if( 2*locM > n )
    locM = n - locM;
  
  if( locM == 0 )
    return 1;
  
  int outVal = n;
  
  int ii;
  //binomailCoeeficient(n, m) = ((((((n * (n-1)) / 2) * (n-2)) / 3) ... ) * (n+1-m)) / m
  //This arrangment gives a whole number for every divide.
  for(ii = 1; ii < locM; ii++)
  {
    outVal = outVal * (n - ii);
    outVal = outVal / (ii+1);
  }
  
  return outVal;
}

int NURBS_T::degreeElevateSurface(std::vector<double>& oldSKnotVec,
                                  const int sDegree,
                                  std::vector<double>& oldTKnotVec,
                                  const int tDegree,
                                  const int Dim,
                                  std::vector<double>& oldCtrlPts,
                                  const int addDegrees,
                                  const char dir)
{
  std::vector<double> newSKnots;
  std::vector<double> newTKnots;
  std::vector<double> newPts;
  
  if( 0 != degreeElevateSurface( oldSKnotVec, sDegree, oldTKnotVec, tDegree,
                                 Dim, oldCtrlPts, addDegrees, dir,
                                 newSKnots, newTKnots, newPts ) )
    return 1;
  
  oldSKnotVec.clear();
  oldTKnotVec.clear();
  oldCtrlPts.clear();
  oldSKnotVec.insert(oldSKnotVec.end(), newSKnots.begin(), newSKnots.end() );
  oldTKnotVec.insert(oldTKnotVec.end(), newTKnots.begin(), newTKnots.end() );
  oldCtrlPts.insert(oldCtrlPts.end(), newPts.begin(), newPts.end() );  

  return 0;
}

int NURBS_T::degreeElevateSurface(const std::vector<double>& oldSKnotVec,
                                  const int sDegree,
                                  const std::vector<double>& oldTKnotVec,
                                  const int tDegree,
                                  const int Dim,
                                  const std::vector<double>& oldCtrlPts,
                                  const int addDegrees,
                                  const char dir,
                                  std::vector<double>& newSKnotVec,
                                  std::vector<double>& newTKnotVec,
                                  std::vector<double>& newCtrlPts)
{
  newSKnotVec.clear();
  newTKnotVec.clear();
  newCtrlPts.clear();
  
  if( 's' == dir )
  {
    //Save unaltered knot vector.
    newTKnotVec.insert( newTKnotVec.end(), oldTKnotVec.begin(), oldTKnotVec.end() );
    
    //Refine each row of points seperately.
    std::vector<double> tempPts;
    std::vector<double> outPts;
    std::vector<double> tempKnotVec;
    int sSize = oldSKnotVec.size() - sDegree - 1;
    
    std::vector<double>::const_iterator cPtIter = oldCtrlPts.begin();
    for( ; cPtIter != oldCtrlPts.end(); cPtIter += (sSize*Dim) )
    {
      tempPts.clear();
      tempPts.insert( tempPts.begin(), cPtIter, (cPtIter+sSize*Dim) );
      
      tempKnotVec.clear();
      outPts.clear();
      if( 0 != degreeElevateCurve( oldSKnotVec, sDegree, Dim, tempPts, addDegrees, 
                                   tempKnotVec, outPts ) )
        return 1;
      
      newCtrlPts.insert( newCtrlPts.end(), outPts.begin(), outPts.end() );
    }
    
    newSKnotVec.insert( newSKnotVec.end(), tempKnotVec.begin(), tempKnotVec.end() );
  }
  else if( 't' == dir )
  {
    //Save unaltered knot vector.
    newSKnotVec.insert( newSKnotVec.begin(), oldSKnotVec.begin(), oldSKnotVec.end() );
    
    int sSize = oldSKnotVec.size() - sDegree - 1;
    int tSize = oldTKnotVec.size() - tDegree - 1;
    
    std::vector<double> tCtrlPts;
    std::vector<double> tempPts;
    std::vector<double> outPts;
    std::vector<double> tempKnotVec;
    
    int ii, jj, dd;
    for( ii = 0; ii < sSize; ii++ )
    {
      tempPts.clear();
      for( jj = 0; jj < tSize; jj++ )
      {
        int offset = (jj*sSize+ii)*Dim;
        for( dd = 0; dd < Dim; dd++ )
          tempPts.push_back( oldCtrlPts[offset+dd] );
      }
      
      tempKnotVec.clear();
      outPts.clear();
      if( 0 != degreeElevateCurve( oldTKnotVec, tDegree, Dim, tempPts, addDegrees,
                                   tempKnotVec, outPts) )
        return 1;
      
      tCtrlPts.insert( tCtrlPts.end(), outPts.begin(), outPts.end() );
    }
    
    newTKnotVec.insert( newTKnotVec.end(), tempKnotVec.begin(), tempKnotVec.end() );
    
    //Transpose the new control points to get them into the correct order.
    tSize = newTKnotVec.size() - (tDegree + addDegrees) - 1;
    for( ii = 0; ii < tSize; ii++ )
    {
      for( jj = 0; jj < sSize; jj++ )
      {
        int offset = (jj*tSize+ii)*Dim;
        std::vector<double>::iterator iter = tCtrlPts.begin() + offset;
        for( dd = 0; dd < Dim; dd++ )
          newCtrlPts.push_back( *(iter+dd) );
      }
    }
  }
  else
  {
    std::cerr << "Incorrect direction parameter given in NURBS_T::degreeElevateSurface()"
    << std::endl;
    return 1;
  }
                
  return 0;
}

int NURBS_T::degreeElevateVolume(std::vector<double>& oldSKnotVec, 
                                 const int sDegree, 
                                 std::vector<double>& oldTKnotVec, 
                                 const int tDegree, 
                                 std::vector<double>& oldUKnotVec, 
                                 const int uDegree, 
                                 const int Dim, 
                                 std::vector<double>& oldCtrlPts, 
                                 const int addDegrees, 
                                 const char dir) 
{ 
  std::vector<double> newSKnots; 
  std::vector<double> newTKnots; 
  std::vector<double> newUKnots; 
  std::vector<double> newPts; 
  
  if( 0 != degreeElevateVolume( oldSKnotVec, sDegree,  
                                oldTKnotVec, tDegree, 
                                oldUKnotVec, uDegree, 
                                Dim, oldCtrlPts, addDegrees, dir, 
                                newSKnots, newTKnots, newUKnots, newPts ) ) 
    return 1; 
  
  oldSKnotVec.clear(); 
  oldTKnotVec.clear(); 
  oldUKnotVec.clear(); 
  oldCtrlPts.clear(); 
  oldSKnotVec.insert(oldSKnotVec.end(), newSKnots.begin(), newSKnots.end() ); 
  oldTKnotVec.insert(oldTKnotVec.end(), newTKnots.begin(), newTKnots.end() ); 
  oldUKnotVec.insert(oldUKnotVec.end(), newUKnots.begin(), newUKnots.end() ); 
  oldCtrlPts.insert(oldCtrlPts.end(), newPts.begin(), newPts.end() );   
  
  return 0; 
} 

int NURBS_T::degreeElevateVolume(const std::vector<double>& oldSKnotVec, 
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
                                 std::vector<double>& newCtrlPts) 
{ 
  newSKnotVec.clear();
  newTKnotVec.clear();
  newUKnotVec.clear();
  newCtrlPts.clear();
  
  std::vector<double> tempInPts;
  std::vector<double> tempOutPts;
  
  //Use surface elavation on layers of control points.
  if( 's' == dir || 't' == dir ) 
  { 
    //Save unchanges knot vector.
    newUKnotVec.clear();
    newUKnotVec.insert( newUKnotVec.begin(), oldUKnotVec.begin(), oldUKnotVec.end() );
    
    int uLayers = oldUKnotVec.size() - uDegree - 1;
    int sCols = oldSKnotVec.size() - sDegree - 1;
    int tRows = oldTKnotVec.size() - tDegree - 1;
    
    int layer;
    for( layer = 0; layer < uLayers; layer++ )
    {
      int row, col;
      for( row = 0; row < tRows; row++ )
      {
        for( col = 0; col < sCols; col++ )
        {
          int ii;
          for( ii = 0; ii < Dim; ii++ )
            tempInPts.push_back( oldCtrlPts[Dim*(layer*sCols*tRows+row*sCols+col)+ii] );
        }
      }
      
      //Coll degree elevation on the surface.  Last loop will set new knot vectors for S and T.
      newSKnotVec.clear(); 
      newTKnotVec.clear();
      if( 0 != degreeElevateSurface(oldSKnotVec, sDegree, oldTKnotVec, tDegree,
                                    Dim, tempInPts, addDegrees, dir, newSKnotVec,
                                    newTKnotVec, tempOutPts) )
        return 1;
      
      newCtrlPts.insert( newCtrlPts.end(), tempOutPts.begin(), tempOutPts.end() );
      
      tempInPts.clear();
      tempOutPts.clear();
    }
  } 
  else if( 'u' == dir ) 
  { 
    //Save S and T knot vectors.
    newSKnotVec.insert( newSKnotVec.end(), oldSKnotVec.begin(), oldSKnotVec.end() );
    newTKnotVec.insert( newTKnotVec.end(), oldTKnotVec.begin(), oldTKnotVec.end() );
    
    std::vector<double> tempKnotVec;
    std::vector<double> tCtrlPts;
    
    int tLayers = oldTKnotVec.size() - tDegree - 1;
    int sCols = oldSKnotVec.size() - sDegree - 1;
    int uRows = oldUKnotVec.size() - uDegree - 1;
    
    int layer;
    for( layer = 0; layer < tLayers; layer++ )
    {
      int row, col;
      for( row = 0; row < uRows; row++ )
      {
        for( col = 0; col < sCols; col++ )
        {
          int ii;
          int offset = (layer*sCols+row*sCols*tLayers+col)*Dim;
          std::vector<double>::const_iterator iter = oldCtrlPts.begin() + offset;
          for( ii = 0; ii < Dim; ii++ )
            tempInPts.push_back( *(iter+ii) );
        }
      }
      
      //Coll degree elevation on the surface.  Last loop will set new knot vectors for  U.
      tempKnotVec.clear();
      newUKnotVec.clear();
      if( 0 != degreeElevateSurface(oldSKnotVec, sDegree, oldUKnotVec, uDegree,
                                    Dim, tempInPts, addDegrees, 't', tempKnotVec,
                                    newUKnotVec, tempOutPts ) )
        return 1;
      
      tCtrlPts.insert( tCtrlPts.end(), tempOutPts.begin(), tempOutPts.end() );
      
      tempInPts.clear();
      tempOutPts.clear();      
    }
    
    //Reorder the control points.
    int sSize = newSKnotVec.size() - sDegree - 1;
    int tSize = newTKnotVec.size() - tDegree - 1;
    int uSize = newUKnotVec.size() - (uDegree + addDegrees) - 1;
    int ii, jj, kk, dd;
    for( ii = 0; ii < uSize; ii++ )
    {
      for( jj = 0; jj < tSize; jj++ )
      {
        for( kk = 0; kk < sSize; kk++ )
        {
          int offset = (ii*sSize+jj*sSize*uSize+kk)*Dim;
          std::vector<double>::iterator iter = tCtrlPts.begin() + offset;
          for( dd = 0; dd < Dim; dd++ )
            newCtrlPts.push_back( *(iter+dd) );
        }
      }
    }
  } 
  else 
  { 
    std::cerr << "Incorrect direction parameter given in NURBS_T::degreeElevateVolume()" 
    << std::endl; 
    return 1; 
  } 
  return 0; 
}

int NURBS_T::curvePoint(const std::vector<double>& knotVec,
                        const int cDegree,
                        const int Dim,
                        const std::vector<double>& ctrlPts,
                        const double s,
                        std::vector<double>& outPts)
{
  std::vector<double> basis;
  
  double span = basisFunctions( knotVec, cDegree, s, basis );
  if( span < 0 )
    return 1;
  
  outPts.clear();
  outPts.resize( Dim, 0.0 );
  
  int offset = (span-cDegree)*Dim;
  int dd, ii;
  for( ii = 0; ii <= cDegree; ii++ )
    for( dd = 0; dd < Dim; dd++ )
      outPts[dd] += basis[ii]*ctrlPts[offset+ii*Dim+dd];
  
  return 0;
}

int NURBS_T::surfacePoint(const std::vector<double>& sKnots,
                          const int sDegree,
                          const std::vector<double>& tKnots,
                          const int tDegree,
                          const int Dim,
                          const std::vector<double>& ctrlPts,
                          const double s,
                          const double t,
                          std::vector<double>& outPts)
{
  std::vector<double> sBasis;
  std::vector<double> tBasis;
  
  double sSpan = basisFunctions( sKnots, sDegree, s, sBasis );
  if( sSpan < 0 )
    return 1;
  
  double tSpan = basisFunctions( tKnots, tDegree, t, tBasis );
  if( tSpan < 0 )
    return 1;
  
  outPts.clear();
  outPts.resize( Dim, 0.0 );
  double* temp = new double[Dim];
  int dd, ii, jj;
  
  int sind = sSpan - sDegree;
  int sPts = sKnots.size() - sDegree - 1;
  
  for( ii = 0; ii <= tDegree; ii++ )
  {
    for( dd = 0; dd < Dim; dd++ )
      temp[dd] = 0.0;
    
    int tind = tSpan - tDegree + ii;
    int offset = (tind*sPts+sind)*Dim;
    
    for( jj = 0; jj <= sDegree; jj++ )
      for( dd = 0; dd < Dim; dd++ )
        temp[dd] += sBasis[jj] * ctrlPts[offset+jj*Dim+dd];
    
    for( dd = 0; dd < Dim; dd++ )
      outPts[dd] += tBasis[ii]*temp[dd];
  }
  
  delete [] temp;
  return 0;
}

//EOF

