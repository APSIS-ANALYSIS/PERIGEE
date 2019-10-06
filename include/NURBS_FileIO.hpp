#ifndef NURBS_FILE_IO_HPP
#define NURBS_FILE_IO_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "vtkVersion.h"

class vtkUnstructuredGrid;
class vtkDoubleArray;
class BezierElem;

//!NUBS Tools namespace
namespace NURBS_T
{
//! Bezier element output data.
/*!
  This class is used to pass element data to the file output
  routines.  Data can be stored at element parametric values
  and can include scalars, vectors and tensors (pretty much
  what the vtk file formats can store).
 */
    class BezElemData
    {
  public:
      BezElemData();
      ~BezElemData();
      
      int initPointDataVector(const char* name, 
                              const int numComponents);
      int setPointDataVectorComponent(const char* name, 
                                      const int gId,
                                      const int component, 
                                      const double value);
      
   /*   
  private:
        friend int NURBS_T::writeVolume(const char* fileName,
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
                                        const bool useXML);
     */ 
        std::map<const std::string, vtkDoubleArray*> mPointVectors;
    };

    
//! Read a NURBS file
/*!
  Reads a NURBS file.  This function does not do a lot of error checking
  so the format of the input file needs to be correct or the results may not
  be valid.
  \param infile The input file.  Should be in the same state as if it was
         just opened.
  \param sKnots Returns the S knot values.
  \param tKnots Returns the T knot values.
  \param uKnots Returns the U knot values.
  \param sDegree Returns the polynomail degree in the S direction.
  \param tDegree Returns the polynomail degree in the T direction.
  \param uDegree Returns the polynomail degree in the U direction.
  \param numCPts Returns the number of control points.
  \param ctrlPts Returns the coordinates of the control points.  Coordinates
         are stored as {x1,y1,z1,w1,x2,y2,z2,w2,...,xn,yn,zn,wn}.
  \return 0 if no errors during input, 1 otherwise.       
 */
    int readFile(std::ifstream& infile,
                 std::vector<double>& sKnots,
                 std::vector<double>& tKnots,
                 std::vector<double>& uKnots,
                 int& sDegree, int& tDegree,
                 int& uDegree, int& numCPts,
                 std::vector<double>& ctrlPts);
  
//! Writes the control net to a vtk file.
/*!
  \param fileName Base name for output file.  Will be appended with .vtk.
  \param sKnots S knot vector
  \param tKnots T knot vector
  \param uKnots U Knot Vector.  Should contain single 0 for surfaces.
  \param sDegree Polynomial degree in S direction
  \param tDegree Polynomial degree in T direction
  \param uDegree Polynomial degree in U direction.  Should be zero for surfaces.
  \param ctrlPts Coordinates of the control points.  Coordinates
         are stored as {x1,y1,z1,w1,x2,y2,z2,w2,...,xn,yn,zn,wn}
         (no z values for 2D).
  \param spatialDim Dimension if writing a surface (2 or 3).
  \return 0 if no errors during output, 1 otherwise.
 */
    int writeControlNet(const char* fileName,
                        const std::vector<double>& sKnots,
                        const std::vector<double>& tKnots,
                        const std::vector<double>& uKnots,
                        const int sDegree,
                        const int tDegree,
                        const int uDegree,
                        const std::vector<double>& ctrlPts,
                        const int spatialDim);

//! Writes the constrol net of a 2D curve to a file.
/*!
  This function writes the control net of a 2D curve to a file.
  Output is (x,y) pairs in column format of the control points.
  \param fileName Base name of file to be writen.  Will be appended with .txt.
  \param ctrlPts Coordinates of the control points.  Coordinates should be
         stored as {x1,y1,w1,...,nx,yn,wn}.
  \return 0 fi no errors, 1 otherwise.       
 */
    int writeControlNet2D(const char* fileName,
                          const std::vector<double>& ctrlPts);


//! Writes a 2D curve to a file.
/*!
  This function writes a 2D curve to a file.
  Output is (x,y) pairs in column format.
  \param fileName Base name of file to be writen.  Will be appended with .txt.
  \param knots The knot vector
  \param cDegree The polynomial degree of the curve.
  \param ctrlPts Coordinates of the control points.  Coordinates should be
         stored as {x1,y1,w1,...,nx,yn,wn}.
  \param numSegments The number of line segements to write per bezier section.
  \return 0 fi no errors, 1 otherwise.       
 */
    int write2DCurve(const char* fileName,
                     const std::vector<double>& knots,
                     const int cDegree,
                     const std::vector<double>& ctrlPts,
                     const int numSegments = 10);

//! Writes a surface to a vtk file.
/*!
  This function writes a NURBS surface to a vtk file.
  \param fileName Base name of file to be writen.  Will be appended with .vtk.
  \param sKnots The knot vector in the S direction.
  \param sDegree The polynomial degree in the S direction.
  \param tKnots The knot vector in the T direction.
  \param tDegree The polynomizl degree  in the T direction.
  \param ctrlPts Coordinates of the homogeneous Bezier control points.  Coordinates should be
         stored as {x1,y1,z1,w1,...,nx,yn,zn,wn}.
  \param dispCtrl The homogeneous displacements of the control points.
         Note, this array should also contain the weights of each point
         (ordered same as ctrlPts).  If no displacements need to be written
         pass in an empty vector.
  \param spatialDim Dimension if writing a surface (2 or 3).         
  \param numSegments The number of devisions in the s and t direction to write
         per bezier section.
  \return 0 fi no errors, 1 otherwise.       
 */    
    int writeSurface(const char* fileName,
                     const std::vector<double>& sKnots,
                     const int sDegree,
                     const std::vector<double>& tKnots,
                     const int tDegree,
                     const std::vector<double>& ctrlPts,
                     const std::vector<double>& dispCtrl,
                     const int spatialDim,
                     const int numSegments = 10);

//! Writes a surface to a vtk file.
/*!
  This function writes a NURBS surface that has been decomposed into
  Bezier segments to a vtk file.
  \param fileName Base name of file to be writen.  Will be appended with .vtk.
  \param bSegments The list of Bezier segments to write.
  \param sKnots The knot vector in the S direction.
  \param sDegree The polynomial degree in the S direction.
  \param tKnots The knot vector in the T direction.
  \param tDegree The polynomizl degree  in the T direction.
  \param ctrlPts Coordinates of the homogeneous Bezier control points.  Coordinates should be
         stored as {x1,y1,z1,w1,...,nx,yn,zn,wn}.
  \param dispCtrl The homogeneous displacements of the control points.
         Note, this array should also contain the weights of each point
         (ordered same as ctrlPts).  If no displacements need to be written
         pass in an empty vector.
  \param spatialDim Dimension if writing a surface (2 or 3).         
  \param numSegments The number of devisions in the s and t direction to write
         per bezier section.
  \return 0 fi no errors, 1 otherwise.       
*/
    int writeSurface(const char* fileName,
                     const std::vector<BezierElem*>& bSegments,
                     const std::vector<double>& sKnots,
                     const int sDegree,
                     const std::vector<double>& tKnots,
                     const int tDegree,
                     const std::vector<double>& ctrlPts,
                     const std::vector<double>& dispCtrl,
                     const int spatialDim,
                     const int numSegments = 10);

//! Writes a NURBS volume to a vtk file.
/*!
  This function writes a NURBS volume that has been decomposed into
  Bezier segments to a vtk file.
  \param fileName Base name of file to be writen.  Will be appended with .vtk
         for legacy format or .vtu for XML format depending on useXML flag.
  \param bSegments The list of Bezier segments to write.
  \param sKnots The knot vector in the S direction.
  \param sDegree The polynomial degree in the S direction.
  \param tKnots The knot vector in the T direction.
  \param tDegree The polynomizl degree  in the T direction.
  \param uKnots The knot vector in the U direction.
  \param uDegree The polynomizl degree in the U direction.
  \param ctrlPts Coordinates of the homogeneous Bezier control points.  Coordinates should be
         stored as {x1,y1,z1,w1,...,nx,yn,zn,wn}.
  \param dispCtrl The homogeneous displacements of the control points.
         Note, this array should also contain the weights of each point
         (ordered same as ctrlPts).  If no displacements need to be written
         pass in an empty vector.
  \param elemData Data that has been calculated on the element.  Element data
         is interpolated to the output element nodes.  See documentation for
         BezElemData class for more information.
  \param numSegmentsS The number of devisions per element in the s direction.
  \param numSegmentsT The number of devisions per element in the t direction.
  \param numSegmentsU The number of devisions per element in the u direction.
  \param useXML Set to true to write XML based vtk files, false to write legacy
         vtk files.
  \return 0 fi no errors, 1 otherwise.       
*/
    int writeVolume(const char* fileName,
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
                    const int numSegmentsS = 1,
                    const int numSegmentsT = 1,
                    const int numSegmentsU = 1,
                    const bool useXML = false);
    
/////////////////////////////////////////////////////////////////
////  Place Helper Functions below this line ////
/////////////////////////////////////////////////////////////////
    
// Helper function for reading input files.    
    int readKnotVector(std::istringstream& sstrm,
                       std::vector<double>& knots);
    
// Sets point data for vtk file.
    int setPoints(const std::vector<double>& ctrlPts,
                  const int spatialDim,
                  vtkUnstructuredGrid* gridData,
                  const bool writeWeight = false);

// Set the displacements for the points.
    int setDisp(const std::vector<double>& dispPts,
                const int spatialDim,
                vtkUnstructuredGrid* gridData);
    
    //Set vector data on the points.
    int setPtVecData(const char* vecName,
                     const std::vector<double> data,
                     const int dim,
                     vtkUnstructuredGrid* gridData);
    
// Sets element data for vtk file.
    int setElements( const int numRows,
                     const int numCols,
                     const int numLayers,
                     const int spatialDim,
                     vtkUnstructuredGrid* gridData);

// Write the vtk file.
    int write(const char* fileName,
              vtkUnstructuredGrid* gridData,
              const bool useXML);    
}
#endif
