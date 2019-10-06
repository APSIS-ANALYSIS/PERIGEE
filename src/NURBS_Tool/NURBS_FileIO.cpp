#include "NURBS_Bezier.hpp"
#include "NURBS_FileIO.hpp"
#include "NURBS_Tools.hpp"

#include <math.h>
#include <assert.h>
#include <set>
#include <iomanip>

#include "vtkDoubleArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkQuad.h"
#include "vtkHexahedron.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"


NURBS_T::BezElemData::BezElemData()
{
}

NURBS_T::BezElemData::~BezElemData()
{
  std::map<const std::string, vtkDoubleArray*>::iterator pvIter;
  pvIter = mPointVectors.begin();
  for( ; pvIter != mPointVectors.end(); pvIter++ )
    pvIter->second->Delete();
}

int NURBS_T::BezElemData::initPointDataVector(const char* name, 
    const int numComponents)
{
  //Make sure not adding duplicates.
  std::map<const std::string, vtkDoubleArray*>::iterator iter;
  iter = mPointVectors.begin();
  for( ; iter != mPointVectors.end(); iter++ )
  {
    if( iter->first.compare( name ) == 0 )
    {
      std::cerr << "ERROR: duplicate name passed into NURBS_T::BezElemData::initPointDataVector()" << std::endl
        << "Aborting initializiation" << std::endl;
      return 1;
    }
  }
  std::string vName(name);

  vtkDoubleArray* dVec = vtkDoubleArray::New();
  dVec->SetNumberOfComponents( numComponents );
  dVec->SetName( name );

  mPointVectors[vName] = dVec;

  return 0;
}

int NURBS_T::BezElemData::setPointDataVectorComponent(const char* name, 
    const int gId,
    const int component, 
    const double value)
{
  //Find the vector.
  std::map<const std::string, vtkDoubleArray*>::iterator iter;
  iter = mPointVectors.begin();
  for( ; iter != mPointVectors.end(); iter++ )
  {
    if( iter->first.compare( name ) == 0 )
      break;
  }

  if( iter == mPointVectors.end() )
  {
    std::cerr << "ERROR: could not find vector (" << name << ") in NURBS_T::BezElemData::addPointDataVectorComponent()" << std::endl
      << "Aborting..." << std::endl;
    return 1;
  }

  //Now make sure the component is a valid index.
  int numComponents = iter->second->GetNumberOfComponents();
  if( component < 0 || component >= numComponents )
  {
    std::cerr << "Invalid compenent value (" << component << ") for vector (" << name 
      << ") in NURBS_T::BezElemData::addPointDataVectorComponent(): " << component << std::endl
                                                                         << "Aborting..." << std::endl;
    return 1;
  }

  iter->second->InsertComponent( gId, component, value );

  return 0;
}

int NURBS_T::readFile(std::ifstream& infile,
    std::vector<double>& sKnots,
    std::vector<double>& tKnots,
    std::vector<double>& uKnots,
    int& sDegree, int& tDegree,
    int& uDegree, int& numCPts,
    std::vector<double>& ctrlPts)
{
  //Get the initial information
  std::istringstream sstrm;
  std::string line;
  std::string sWord;
  char inWord[1001];
  char x;

  int hits = 0;
  while( hits < 8 && !infile.eof() )
  {
    infile.getline(inWord, 1000);
    line = inWord;
    sstrm.str( line );
    sstrm >> inWord;

    sWord = inWord;
    if( sWord == "TYPE" )
    {
      sstrm >> x;
      if( x != '=' )
      {
        std::cerr << "Bad file format: 1" << std::endl;
        return 1;
      }
      sstrm >> inWord;
      sWord = inWord;
      if( sWord != "NURBS" )
      {
        std::cerr << "The file format is not NURBS.  Cannot read file." << std::endl;
        return 1;
      }

      hits++;
    }
    else if( sWord == "GLOBAL_S" )
    {
      readKnotVector(sstrm, sKnots);
      hits++;
    }
    else if( sWord == "GLOBAL_T" )
    {
      readKnotVector(sstrm, tKnots);
      hits++;
    }    
    else if( sWord == "GLOBAL_U" )
    {
      readKnotVector(sstrm, uKnots);
      hits++;
    }
    else if( sWord == "DEGREE_S" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.1" << std::endl;
        return 1;
      }
      sstrm >> sDegree;
      hits++;
    }
    else if( sWord == "DEGREE_T" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.2" << std::endl;
        return 1;
      }
      sstrm >> tDegree;
      hits++;
    }    
    else if( sWord == "DEGREE_U" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.3" << std::endl;
        return 1;
      }
      sstrm >> uDegree;
      hits++;
    }
    else if( sWord == "NUM_CP" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 4.1" << std::endl;
        return 1;
      }
      sstrm >> numCPts;
      hits++;
    }

    sstrm.clear();
  }

  if( 8 != hits )
  {
    std::cerr << "Error reading file.  Something is wrong with the format\n" << std::endl;
    return 1;
  }

  std::vector<double>::iterator iter;
  std::cout << "S Knots: ";
  for( iter = sKnots.begin(); iter != sKnots.end(); iter++ )
    std::cout << *iter << " ";
  std::cout << std::endl << "T Knots: ";
  for( iter = tKnots.begin(); iter != tKnots.end(); iter++ )
    std::cout << *iter << " ";
  std::cout << std::endl << "U Knots: ";
  for( iter = uKnots.begin(); iter != uKnots.end(); iter++ )
    std::cout << *iter << " ";
  std::cout << std::endl;

  std::cout << "S Degree: " << sDegree << std::endl;
  std::cout << "T Degree: " << tDegree << std::endl;
  std::cout << "U Degree: " << uDegree << std::endl;
  std::cout << "Num C Pts: " << numCPts << std::endl;

  //Read in the control points.
  int ii;
  for( ii = 0; ii < numCPts*4; ii++ )
  {
    if( infile.eof() )
    {
      std::cerr << "Wrong number of control points." << std::endl;
      return 1;
    }

    double coord;
    infile >> coord;
    ctrlPts.push_back( coord );
  }

  if( (int)ctrlPts.size() != numCPts*4 )
  {
    std::cerr << "Error reading control points" << std::endl
      << ctrlPts.size() << std::endl;
    return 1;
  }

  return 0;
}

int NURBS_T::readKnotVector(std::istringstream& sstrm,
    std::vector<double>& knots)
{
  char x;

  //Discard = and [ characters
  sstrm >> x;
  if( x != '=' )
  {
    std::cerr << "Bad file format: 2.1" << std::endl;
    return 1;
  }
  sstrm >> x;
  if( x != '[' )
  {
    std::cerr << "Bad file format: 2.2" << std::endl;
    return 1;
  }

  sstrm >> x;
  while( x != ']' && !sstrm.eof() )
  {
    sstrm.putback( x );
    double knot;
    sstrm >> knot;
    knots.push_back( knot );
    sstrm >> x;
  }

  return 0;
}

int NURBS_T::writeControlNet(const char* fileName,
    const std::vector<double>& sKnots,
    const std::vector<double>& tKnots,
    const std::vector<double>& uKnots,
    const int sDegree,
    const int tDegree,
    const int uDegree,
    const std::vector<double>& ctrlPts,
    const int spatialDim)
{
  vtkUnstructuredGrid* gridData = vtkUnstructuredGrid::New();

  if( 0 != setPoints( ctrlPts, spatialDim, gridData, true ) )
    return 1;

  const int numRows = tKnots.size() - tDegree - 1;
  const int numCols = sKnots.size() - sDegree - 1;
  const int numLayers = 0 == uDegree ? 1 : uKnots.size() - uDegree - 1;

  if( 0 != setElements( numRows, numCols, numLayers, spatialDim, gridData ) )
    return 1;

  if( 0 != write( fileName, gridData, false ) )
    return 1;

  gridData->Delete();

  return 0;
}

int NURBS_T::setPoints(const std::vector<double>& ctrlPts,
    const int spatialDim,
    vtkUnstructuredGrid* gridData,
    const bool writeWeight)
{
  const int cpDim = spatialDim + 1; //Account for weight.
  const int numPts = ctrlPts.size() / cpDim;
  double xyz[3] = {0,0,0};

  vtkPoints* points = vtkPoints::New();

  vtkDoubleArray* dVector = NULL;
  if( true == writeWeight )
  {
    dVector = vtkDoubleArray::New();
    dVector->SetNumberOfComponents(1);
    dVector->SetName("Weight");
  }

  int ii;
  for( ii = 0; ii < numPts; ii++ )
  {
    switch( spatialDim )
    {
      case 3:
        xyz[2] = ctrlPts[ii*cpDim+2];
      case 2:
        xyz[1] = ctrlPts[ii*cpDim+1];
      default:
        break;
    }
    xyz[0] = ctrlPts[ii*cpDim];

    points->InsertPoint( ii, xyz );
    if( true == writeWeight )
      dVector->InsertValue( ii, ctrlPts[ii*cpDim + spatialDim] ); 
  }

  gridData->SetPoints( points );
  points->Delete();

  if( true == writeWeight )
  {
    gridData->GetPointData()->AddArray( dVector );
    dVector->Delete();
  }
  return 0;
}

int NURBS_T::setElements( const int numRows,
    const int numCols,
    const int numLayers,
    const int spatialDim,
    vtkUnstructuredGrid* gridData)
{
  if( 1 == numLayers )
  {
    //NURBS surface
    int row, col;
    for( row = 0; row < numRows-1; row++ )
    {
      for( col = 0; col < numCols-1; col++ )
      {
        vtkCell* cell = vtkQuad::New();
        cell->GetPointIds()->SetId( 0, row*numCols + col);
        cell->GetPointIds()->SetId( 1, row*numCols + col+1);
        cell->GetPointIds()->SetId( 2, (row+1)*numCols + col+1);
        cell->GetPointIds()->SetId( 3, (row+1)*numCols + col);

        gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
        cell->Delete();
      }
    }
  }
  else
  {
    //Trivariate NURBS
    int row, col, layer;
    for( layer = 0; layer < numLayers-1; layer++ )
    {
      for( row = 0; row < numRows-1; row++ )
      {
        for( col = 0; col < numCols-1; col++ )
        {
          vtkCell* cell = vtkHexahedron::New();

          cell->GetPointIds()->SetId( 0, layer*numRows*numCols+row*numCols + col);
          cell->GetPointIds()->SetId( 1, layer*numRows*numCols+row*numCols + col+1);
          cell->GetPointIds()->SetId( 2, layer*numRows*numCols+(row+1)*numCols + col+1);
          cell->GetPointIds()->SetId( 3, layer*numRows*numCols+(row+1)*numCols + col);

          cell->GetPointIds()->SetId( 4, (layer+1)*numRows*numCols+row*numCols + col);
          cell->GetPointIds()->SetId( 5, (layer+1)*numRows*numCols+row*numCols + col+1);
          cell->GetPointIds()->SetId( 6, (layer+1)*numRows*numCols+(row+1)*numCols + col+1);
          cell->GetPointIds()->SetId( 7, (layer+1)*numRows*numCols+(row+1)*numCols + col);

          gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
          cell->Delete();
        }
      }
    }
  }

  return 0;
}

int NURBS_T::write(const char* fileName,
    vtkUnstructuredGrid* gridData,
    const bool useXML)
{
  std::string fName( fileName );
  if( useXML )
  {
    fName.append( ".vtu" );

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName( fName.c_str() );
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput( gridData );
#else
    writer->SetInputData( gridData );
#endif
    if( 1 != writer->Write() )
    {
      std::cerr << "ERROR: Failed to write results file." << std::endl;
      writer->Delete();
      return 1;
    }  

    writer->Delete();    

  }
  else
  {
    fName.append( ".vtk" );

    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    writer->SetFileName( fName.c_str() );
#if VTK_MAJOR_VERSION <= 5    
    writer->SetInput( gridData );
#else
    writer->SetInputData( gridData );
#endif    
    if( 1 != writer->Write() )
    {
      std::cerr << "ERROR: Failed to write results file." << std::endl;
      writer->Delete();
      return 1;
    }  

    writer->Delete(); 
  }
  return 0;
}

int NURBS_T::writeControlNet2D(const char* fileName,
    const std::vector<double>& ctrlPts)
{
  std::string fName( fileName );
  fName.append( ".txt" );

  //Open the file
  std::ofstream outfile( fName.c_str(), std::ifstream::out);
  if( !outfile )
  {
    std::cerr << "Unable to open file to write 2D control net: " << fName << std::endl;
    return 1;
  }

  int ii;
  for( ii = 0; ii < (int)ctrlPts.size() / 3; ii++ )
    outfile << ctrlPts[3*ii] << " " << ctrlPts[3*ii+1] << std::endl;

  outfile.close();

  return 0;
}

int NURBS_T::write2DCurve(const char* fileName,
    const std::vector<double>& knots,
    const int cDegree,
    const std::vector<double>& ctrlPts,
    const int numSegments)
{
  std::string fName( fileName );
  fName.append( ".txt" );

  //Open the file
  std::ofstream outfile( fName.c_str(), std::ifstream::out);
  if( !outfile )
  {
    std::cerr << "Unable to open file to write 2D curve: " << fName << std::endl;
    return 1;
  }

  //Get the bezier decomposition for the curve.
  std::vector<double> bCtrlPts;
  std::vector<NURBS_T::BezierElem*> bSegments;

  if( 0 != NURBS_T::decomposeCurve( knots, cDegree, ctrlPts, 3,
        bCtrlPts, bSegments ) )
  {
    std::cerr << "ERROR" << std::endl;
    return 1;
  }

  std::vector<double> outPoints;

  std::vector<double>::iterator iter;
  iter = bCtrlPts.begin();
  double step = 1.0 / numSegments;
  int ii, jj, mm;
  for( jj = 0; jj < (int)bSegments.size(); jj++ )
  {
    //Get the control points for this Bezier.
    std::vector<double> bPts;
    for( ii = 0; ii < cDegree+1; ii++ )
    {
      for( mm = 0; mm < 3; mm++ )
      {
        bPts.push_back(*iter);
        iter++;
      }
    }
    iter -= 3;

    //Get the new points.
    std::vector<double> pt;
    double u = step;
    for( ii = 1; ii <= numSegments; ii++ )
    {
      if( 0 != deCasteljauCurve( bPts, 3, u, pt ) )
        return 1;

      outPoints.insert( outPoints.end(), pt.begin(), pt.end() );

      u += step;
    }
  }

  projectDown( outPoints, 2 );

  //Write the first control point.  It is not contained in outPoints.
  outfile << ctrlPts[3*bSegments[0]->gPtsMap[0]]
    << " " << ctrlPts[3*bSegments[0]->gPtsMap[0] + 1] << std::endl;

  //Write the points to file.
  int size = outPoints.size() / 3;
  for( ii = 0; ii < size; ii++ )
    outfile << outPoints[3*ii] << " " << std::setprecision(12) << outPoints[3*ii+1] << std::endl;

  //Delete the Bezier segments.
  while( !bSegments.empty() )
  {
    delete bSegments.back();
    bSegments.pop_back();
  }

  return 0;
}

int NURBS_T::writeSurface(const char* fileName,
    const std::vector<BezierElem*>& bSegments,
    const std::vector<double>& sKnots,
    const int sDegree,
    const std::vector<double>& tKnots,
    const int tDegree,
    const std::vector<double>& ctrlPts,
    const std::vector<double>& dispCtrl,
    const int spatialDim,
    const int numSegments)
{
  const int cpDim = spatialDim+1;

  //Get the number of segments in each direction.
  //(ie. number of unique knot values minus 1).
  int ii;
  int numSKnots = 0;
  int numTKnots = 0;

  for( ii = 0; ii < (int)sKnots.size()-1; ii++ )
  {
    numSKnots++;
    while( ii < (int)sKnots.size()-1 && sKnots[ii] == sKnots[ii+1] )
      ii++;
  }
  for( ii = 0; ii < (int)tKnots.size()-1; ii++ )
  {
    numTKnots++;
    while( ii < (int)tKnots.size()-1 && tKnots[ii] == tKnots[ii+1] )
      ii++;
  }

  const int numSegsS = numSKnots - 1;
  const int numSegsT = numTKnots - 1;

  assert( (int)bSegments.size() == numSegsS*numSegsT );

  const bool includeDisp = dispCtrl.size() > 0 ? true : false;
  if( includeDisp && dispCtrl.size() != ctrlPts.size() )
  {
    std::cerr << "ERROR (in NURBS_FileIO::writeSurface()): point displacements do not match the "
                                                           << "number of control points" << std::endl
                                                           << "Output file not written." << std::endl;
    return 1;
  }

  std::vector<double> outPoints;
  std::vector<double> outDisp;
  const int numSPts = numSegsS*numSegments+1;
  const int numTPts = numSegsT*numSegments+1;

  //Resize outPoints so that it can be indexed.
  outPoints.resize(cpDim*numSPts*numTPts, 0.0 );
  if( includeDisp )
    outDisp.resize( outPoints.size(), 0.0 );


  int sIndex, tIndex;
  const double step = 1.0 / numSegments;
  const int numBElemPts = (sDegree+1)*(tDegree+1);
  for( tIndex = 0; tIndex < numSegsT; tIndex++ )
  {
    for( sIndex = 0; sIndex < numSegsS; sIndex++ )
    {
      BezierElem* elem = bSegments[tIndex*numSegsS+sIndex];

      std::vector<double> bezPts;
      std::vector<double> dPts;

      //Get the Bezier points.
      for( ii = 0; ii < numBElemPts; ii++ )
      {
        int jj;
        for( jj = 0; jj < cpDim; jj++ )
        {
          bezPts.push_back( ctrlPts[ cpDim*elem->lPtsMap[ii]+jj ] );

          if( includeDisp )
            dPts.push_back( dispCtrl[ cpDim*elem->lPtsMap[ii]+jj ] );
        }
      }


      std::vector<double> pts;
      std::vector<double> ptsD;

      int startS, startT;
      int ss, tt;
      double u, v;
      if( 0 == sIndex )
      {
        startS = 0;
        u = 0;
      }
      else
      {
        startS = 1;
        u = step;
      }      
      if( 0 == tIndex )
      {
        startT = 0;
        v = 0;
      }
      else
      {
        startT = 1;
        v = step;
      }

      std::vector<double> sVec;
      std::vector<double> tVec;

      for( ss = startS; ss <= numSegments; ss++ )
      {
        sVec.push_back( u );
        u += step;
      }
      for( tt = startT; tt <= numSegments; tt++ )
      {
        tVec.push_back( v );
        v += step;
      }

      if( 0 != deCasteljauSurface( bezPts, sDegree, tDegree, cpDim,
            sVec, tVec, pts ) )
        return 1;

      if( includeDisp )
      {
        if( 0 != deCasteljauSurface( dPts, sDegree, tDegree, cpDim, 
              sVec, tVec, ptsD ) )
          return 1;
      }

      for( tt = startT; tt <= numSegments; tt++ )
      {
        int tOffset = (tIndex*numSegments+tt)*numSPts;
        for( ss = startS; ss <= numSegments; ss++ )
        {
          int offset = cpDim*(tOffset+sIndex*numSegments+ss);
          int ptOffset = cpDim*((tt-startT)*sVec.size()+ss-startS);
          for( ii = 0; ii < cpDim; ii++ )
            outPoints[offset+ii] = pts[ptOffset+ii];

          if( includeDisp )
          {
            for( ii = 0; ii < cpDim; ii++ )
              outDisp[offset+ii] = ptsD[ptOffset+ii];
          }
        }
      }
    }
  }

  projectDown( outPoints, spatialDim );


  vtkUnstructuredGrid* gridData = vtkUnstructuredGrid::New();
  if( 0 != setPoints( outPoints, spatialDim, gridData ) )
    return 1;
  if( includeDisp )
  {
    projectDown( outDisp, spatialDim );
    if( 0 != setDisp( outDisp, spatialDim, gridData ) )
      return 1;
  }
  if( 0 != setElements( numTPts, numSPts, 1, spatialDim, gridData ) )
    return 1;

  if( 0 != write( fileName, gridData, false ) )
    return 1;

  gridData->Delete();

  return 0;
}

int NURBS_T::writeSurface(const char* fileName,
    const std::vector<double>& sKnots,
    const int sDegree,
    const std::vector<double>& tKnots,
    const int tDegree,
    const std::vector<double>& ctrlPts,
    const std::vector<double>& dispCtrl,
    const int spatialDim,
    const int numSegments)
{
  std::vector<BezierElem*> bSegments;
  std::vector<double> bCtrlPts;

  if( 0 != decomposeSurface( sKnots, sDegree, tKnots, tDegree,
        ctrlPts, spatialDim+1, bCtrlPts, bSegments ) )
    return 1;

  std::vector<double> outPts;
  if( dispCtrl.size() > 0 )
  {
    if( 0 != transformPoints( dispCtrl, bSegments, spatialDim+1, outPts ) )
      return 1;
  }

  if ( 0 != writeSurface( fileName, bSegments, sKnots, sDegree, tKnots, tDegree,
        bCtrlPts, outPts, spatialDim, numSegments ) )
    return 1;

  // valgrind does not recongnize this way to delete the dynamic pointer 
  // array bSegments. I use iterator instead to avoid leak report in
  // valgrind
  // made by Ju Liu.
  //
  //for(int i=0; i<bSegments.size(); i++)
  //  bSegments[i]->~BezierElem();
  std::vector<NURBS_T::BezierElem*>::iterator it;
  for( it = bSegments.begin(); it<bSegments.end(); it++ )
    delete *it;


  return 0;
}

int NURBS_T::writeVolume(const char* fileName,
    const std::vector<BezierElem*>& bSegments,
    const std::vector<double>& sKnots,
    const int sDegree,
    const std::vector<double>& tKnots,
    const int tDegree,
    const std::vector<double>& uKnots,
    const int uDegree,
    const std::vector<double>& ctrlPts,
    const std::vector<double>& dispCtrl,
    const BezElemData& elemData,
    const int numSegmentsS,
    const int numSegmentsT,
    const int numSegmentsU,
    const bool useXML)
{
  const int cpDim = 4;

  //Get the number of segments in each direction.
  //(ie. number of unique knot values minus 1).
  int ii;
  int numSKnots = 0;
  int numTKnots = 0;
  int numUKnots = 0;

  for( ii = 0; ii < (int)sKnots.size()-1; ii++ )
  {
    numSKnots++;
    while( ii < (int)sKnots.size()-1 && sKnots[ii] == sKnots[ii+1] )
      ii++;
  }
  for( ii = 0; ii < (int)tKnots.size()-1; ii++ )
  {
    numTKnots++;
    while( ii < (int)tKnots.size()-1 && tKnots[ii] == tKnots[ii+1] )
      ii++;
  }
  for( ii = 0; ii < (int)uKnots.size()-1; ii++ )
  {
    numUKnots++;
    while( ii < (int)uKnots.size()-1 && uKnots[ii] == uKnots[ii+1] )
      ii++;
  }

  const int numSegsS = numSKnots - 1;
  const int numSegsT = numTKnots - 1;
  const int numSegsU = numUKnots - 1;

  assert( (int)bSegments.size() == numSegsS*numSegsT*numSegsU );

  std::vector<double> outPoints;
  std::vector<double> outDisp;
  const bool includeDisp = dispCtrl.size() > 0 ? true : false;
  if( includeDisp && dispCtrl.size() != ctrlPts.size() )
  {
    std::cerr << "ERROR: point displacements do not match the number of control points" << std::endl
      << "Output file not written." << std::endl;
    return 1;
  }

  const int numSPts = numSegsS*numSegmentsS+1;
  const int numTPts = numSegsT*numSegmentsT+1;
  const int numUPts = numSegsU*numSegmentsU+1;
  const int numTotalPoints = numSPts*numTPts*numUPts;

  //Resize outPoints so that it can be indexed.
  outPoints.resize(cpDim*numTotalPoints, 0.0 );
  if( includeDisp )
    outDisp.resize(outPoints.size(), 0.0 );

  //Create vectors for elemData.
  std::vector< std::vector<double>* > elemPtVecs;  //Global output data
  std::vector< std::vector<double>* > vPts;        //Localalized element data
  std::vector< std::vector<double>* > vOutPts;     //Local output data.

  std::vector<int> elemPtVecsDim;
  std::map<const std::string, vtkDoubleArray*>::const_iterator pvIter;
  pvIter = elemData.mPointVectors.begin();
  for( ; pvIter != elemData.mPointVectors.end(); pvIter++ )
  {
    std::vector<double>* pVec = new std::vector<double>;

    //Resize the vector.  Add one component for the weight.
    int compSize = pvIter->second->GetNumberOfComponents() + 1;
    pVec->resize( numTotalPoints*compSize, 0.0 );
    elemPtVecs.push_back( pVec );
    elemPtVecsDim.push_back( compSize );

    pVec = new std::vector<double>;
    vPts.push_back( pVec );

    pVec = new std::vector<double>;
    vOutPts.push_back( pVec );
  }

  int sIndex, tIndex, uIndex;
  const double stepS = 1.0 / numSegmentsS;
  const double stepT = 1.0 / numSegmentsT;
  const double stepU = 1.0 / numSegmentsU;
  const int numBElemPts = (sDegree+1)*(tDegree+1)*(uDegree+1);
  std::vector<double>::const_iterator c_iter, d_iter;
  std::vector<double> sVec, tVec, uVec;

  for( uIndex = 0; uIndex < numSegsU; uIndex++ )
  {
    for( tIndex = 0; tIndex < numSegsT; tIndex++ )
    {
      for( sIndex = 0; sIndex < numSegsS; sIndex++ )
      {
        BezierElem* elem = bSegments[uIndex*numSegsS*numSegsT+tIndex*numSegsS+sIndex];

        //Get the parametric values
        double vS = stepS;
        double vT = stepT;
        double vU = stepU;
        if( 0 == sIndex )
          sVec.push_back( 0 );
        if( 0 == tIndex )
          tVec.push_back( 0 );
        if( 0 == uIndex )
          uVec.push_back( 0 );

        for( ii = 1; ii <= numSegmentsS; ii++ )
        {
          sVec.push_back( vS );
          vS += stepS;
        }
        for( ii = 1; ii <= numSegmentsT; ii++ )
        {
          tVec.push_back( vT );
          vT += stepT;
        }
        for( ii = 1; ii <= numSegmentsU; ii++ )
        {
          uVec.push_back( vU );
          vU += stepU;
        }

        std::vector<double> bezPts;
        std::vector<double> dPts;

        //Get the Bezier points.
        for( ii = 0; ii < numBElemPts; ii++ )
        {
          c_iter = ctrlPts.begin() + cpDim*elem->lPtsMap[ii];

          if( includeDisp )
            d_iter = dispCtrl.begin() + cpDim*elem->lPtsMap[ii];

          int jj;
          for( jj = 0; jj < cpDim; jj++ )
          {
            bezPts.push_back( *(c_iter+jj) );

            if( includeDisp )
              dPts.push_back( *(d_iter+jj) );
          }

          jj = 0;
          pvIter = elemData.mPointVectors.begin();
          for( ; pvIter != elemData.mPointVectors.end(); pvIter++ )
          {
            double* tuple = new double[ elemPtVecsDim[jj] - 1 ];

            pvIter->second->GetTuple( elem->lPtsMap[ii], tuple );

            vPts[jj]->insert( vPts[jj]->end(), &tuple[0], &tuple[elemPtVecsDim[jj]-1] );

            //Add the weight for this point
            vPts[jj]->push_back( bezPts.back() );            
            delete [] tuple;
            jj++;
          }
        }

        //Get the points.
        std::vector<double> pts;
        std::vector<double> ptsD;
        if( 0 != deCasteljauVolume( bezPts, sDegree, tDegree, uDegree,
              sVec, tVec, uVec, 4, pts ) )
          return 1;

        if( includeDisp )
        {
          if( 0 != deCasteljauVolume( dPts, sDegree, tDegree, uDegree,
                sVec, tVec, uVec, 4, ptsD ) )
            return 1;
        }


        for( ii = 0; ii < (int)vPts.size(); ii++ )
        {
          vOutPts[ii]->clear();
          if( 0 != deCasteljauVolume( *(vPts[ii]), sDegree, tDegree, uDegree,
                sVec, tVec, uVec, elemPtVecsDim[ii], *(vOutPts[ii]) ) )
            return 1;
        }

        int ss, tt, uu;
        int startS = (sIndex == 0) ? 0 : 1;
        int startT = (tIndex == 0) ? 0 : 1;
        int startU = (uIndex == 0) ? 0 : 1;
        for( uu = startU; uu <= numSegmentsU; uu++ )
        {
          int uOffset = (uIndex*numSegmentsU+uu)*numSPts*numTPts;
          for( tt = startT; tt <= numSegmentsT; tt++ )
          {
            int tOffset = (tIndex*numSegmentsT+tt)*numSPts;
            for( ss = startS; ss <= numSegmentsS; ss++ )
            {
              int offset = cpDim*(uOffset+tOffset+sIndex*numSegmentsS+ss);
              int ptOffset = cpDim*((uu-startU)*sVec.size()*tVec.size()
                  +(tt-startT)*sVec.size()+ss-startS);

              for( ii = 0; ii < cpDim; ii++ )
                outPoints[offset+ii] = pts[ptOffset+ii];

              if( includeDisp )
              {
                for( ii = 0; ii < cpDim; ii++ )
                  outDisp[offset+ii] = ptsD[ptOffset+ii];
              }

              for( ii = 0; ii < (int)vOutPts.size(); ii++ )
              {
                offset = elemPtVecsDim[ii]*(uOffset+tOffset+sIndex*numSegmentsS+ss);
                ptOffset = elemPtVecsDim[ii]*((uu-startU)*sVec.size()*tVec.size()
                    +(tt-startT)*sVec.size()+ss-startS);

                std::vector<double>* globalVec = elemPtVecs[ii];
                std::vector<double>* localVec = vOutPts[ii];
                int jj;
                for( jj = 0; jj < elemPtVecsDim[ii]; jj++ )
                  (*globalVec)[offset+jj] = (*localVec)[ptOffset+jj];
              }
            }
          }
        }

        for( ii = 0; ii < (int)vPts.size(); ii++ )
          vPts[ii]->clear();

        sVec.clear();
        tVec.clear();
        uVec.clear();
      }
    }
  }

  projectDown( outPoints, 3 );

  vtkUnstructuredGrid* gridData = vtkUnstructuredGrid::New();
  if( 0 != setPoints( outPoints, 3, gridData ) )
    return 1;
  if( includeDisp )
  {
    projectDown( outDisp, 3 );
    if( 0 != setDisp( outDisp, 3, gridData ) )
      return 1;
  }

  ii = 0;
  pvIter = elemData.mPointVectors.begin();
  for( ; pvIter != elemData.mPointVectors.end(); pvIter++ )
  {
    projectDown( *(elemPtVecs[ii]), elemPtVecsDim[ii] - 1 );

    if( 0 != setPtVecData( pvIter->first.c_str(), *(elemPtVecs[ii]), elemPtVecsDim[ii] - 1, gridData ) )
      return 1;

    ii++;
  }

  if( 0 != setElements( numTPts, numSPts, numUPts, 3, gridData ) )
    return 1;
  if( 0 != write( fileName, gridData, useXML ) )
    return 1;

  //Cleanup Memory
  gridData->Delete();

  while( !elemPtVecs.empty() )
  {
    delete elemPtVecs.back();
    elemPtVecs.pop_back();

    delete vPts.back();
    vPts.pop_back();

    delete vOutPts.back();
    vOutPts.pop_back();
  }

  return 0;
}

int NURBS_T::setDisp(const std::vector<double>& dispPts,
    const int spatialDim,
    vtkUnstructuredGrid* gridData)
{
  const int cpDim = spatialDim + 1;

  const int numPts = dispPts.size() / cpDim;

  vtkDoubleArray* dVector = vtkDoubleArray::New();
  dVector->SetNumberOfComponents(spatialDim);
  dVector->SetName( "Displacement" );

  int ii;
  std::vector<double>::const_iterator iter = dispPts.begin();
  for( ii = 0; ii < numPts; ii++ )
  {
    if( 2 == spatialDim )
    {
      dVector->InsertTuple2( ii, *iter, *(iter+1) );
      iter += cpDim;
    }
    else
    {
      dVector->InsertTuple3( ii, *iter, *(iter+1), *(iter+2) );
      iter += cpDim;
    }
  }

  gridData->GetPointData()->AddArray( dVector );

  dVector->Delete();

  return 0;
}

int NURBS_T::setPtVecData(const char* vecName,
    const std::vector<double> data,
    const int dim,
    vtkUnstructuredGrid* gridData)
{
  const int cpDim = dim + 1;

  const int numPts = data.size() / cpDim;

  vtkDoubleArray* dVector = vtkDoubleArray::New();
  dVector->SetNumberOfComponents( dim );
  dVector->SetName( vecName );

  int ii, jj;
  double* tuple = new double[dim];
  std::vector<double>::const_iterator iter = data.begin();
  for( ii = 0; ii < numPts; ii++ )
  {
    for( jj = 0; jj < dim; jj++ )
      tuple[jj] = *(iter + jj );
    dVector->InsertTuple( ii, tuple );
    iter += cpDim;
  }

  gridData->GetPointData()->AddArray( dVector );

  dVector->Delete();
  delete [] tuple;

  return 0;
}
//EOF

