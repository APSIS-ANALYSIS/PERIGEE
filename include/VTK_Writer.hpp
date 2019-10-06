#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP
// ==================================================================
// VTK_Writer.hpp
//
// This is the class that contains the functions we use to visualize
// solutions as vtk/vtu files.
//
// Date: Dec 17 2013
// ==================================================================
#include <fstream>

#include "IAGlobal_Mesh_Info.hpp"
#include "IAExtractor.hpp"
#include "ALocal_Elem.hpp"
#include "FEAElement_NURBS_3D_der0_v3.hpp"
#include "FEAElement_NURBS_2D_der0.hpp"
#include "IVisDataPrep.hpp"

#include "hdf5.h"

#include "vtkVersion.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkQuad.h"
#include "vtkHexahedron.h"
#include "vtkCellData.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

class VTK_Writer
{
  public:
    VTK_Writer( const IAGlobal_Mesh_Info * const &gmesh_ptr );
    virtual ~VTK_Writer();

    // --------------------------------------------------------------
    // writeOutput:
    // This is the VTK_Writer class's major function that build vtk
    // data structure and write data in vtk format.
    //
    // The writeOutput function is loaded by different parameter input
    // set. For 3D output, there are three Bernstein polynomial input;
    // For 2D output, there are two Bernstein polynomial input.
    //
    // In the writeOutput function, 
    // 1. a vtkUnstructuredGrid object -- gridData is created.
    // 2. build2/3doutput function is called to pass data to gridData.
    //    gridData carries the full information of the interpolated
    //    data for visualization.
    // 3. If necessary, additional data is written for gridData out of
    //    the buildoutput level.
    // 4. writeGirdObject is called to write each subdomain/processor's
    //    .vtu file.
    // 5. writepvtuFile is called to udpate the .pvtu file.
    // 6. gridData is deleted.
    // 6. writepvdFile is called to update the .pvd file.
    // --------------------------------------------------------------
    // writeOutput for 3D visualization. 
    // --------------------------------------------------------------
    void writeOutput( const std::string &epart_file,
        const IAGlobal_Mesh_Info * const &gmesh_ptr,
        const FEANode * const &fnode_ptr,
        const IALocal_meshSize * const &lmsize_ptr,
        const IAExtractor * const &ext_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        const BernsteinBasis_Array * const &bern_x,
        const BernsteinBasis_Array * const &bern_y,
        const BernsteinBasis_Array * const &bern_z,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML ) const;


    // --------------------------------------------------------------
    // writeOutput for 2D visualization. 
    // --------------------------------------------------------------
    void writeOutput( const std::string &epart_file,
        const IAGlobal_Mesh_Info * const &gmesh_ptr,
        const FEANode * const &fnode_ptr,
        const IALocal_meshSize * const &lmsize_ptr,
        const IAExtractor * const &ext_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        const BernsteinBasis_Array * const &bern_x,
        const BernsteinBasis_Array * const &bern_y,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML ) const;



  private:
    // ! R: the array allocated for access basis function values with size
    //      nLocBas
    double * R;

    // ! nLocBas: number of local basis functions in a single element
    int nLocBas;
     
    // --------------------------------------------------------------
    // ! writeGridObject: based on the info provided by gridData (
    //   generated from buildxdoutput), this writes vtu/vtk files on disk
    //   The generated vtk/vtu files will have a name as
    //   baseName_pxxxx.vtk(vtu), xxxx is the cpu number
    void writeGridObject( vtkUnstructuredGrid * gridData,
        const std::string &baseName, const int &rank, const int &size,
        const bool &isXML ) const;

    // --------------------------------------------------------------
    // ! writepvtuFile: write the pvtu file that organize all vtu files
    //   The generated pvtu file will have name:
    //   baseName.pvtu
    void writepvtuFile( vtkUnstructuredGrid * gridData,
        const std::string &baseName, const int &rank, const int &size,
        const bool &isXML) const;
   
    
    // --------------------------------------------------------------
    // ! writepvdFile_Init: write the top lines of pvd file
    void writepvdFile_Init( const std::string &pvdFName ) const;

    // --------------------------------------------------------------
    // ! writepvdFile: append pvd info in each time step writting
    //   the pvd file will have a name as baseName.pvd
    // \para pre_pvtuname: pre_pvtuname.pvtu is the corresponding pvtu
    //                     file
    void writepvdFile( const std::string &baseName, 
       const std::string &pre_pvtuname, const int &rank, const int &size,
       const double &sol_time, const bool &isXML ) const; 

    // --------------------------------------------------------------
    // ! build3doutput: assigns to gridData the sampling points, 
    //                  sampling elements, and data that associated with
    //                  them.
    // \para gridData: grid data holder that stores all info to be written
    // \para epart_file: epart.h5
    // \para nElem: the number of global element
    // \para fNode_ptr: pointer to control points
    // \para lmsize_ptr: local mesh size pointer
    // \para ext_ptr: extraction operators
    // \para lien_ptr: pointer to LIEN
    // \para lelem_ptr: local element index
    // \para vdata_ptr: visualization data struct 
    // \para bern_x _y _z: bernstein polynomials for generating element
    // \para pointArrays: double pointer holding the solution to be vis
    // --------------------------------------------------------------
    void build3doutput(vtkUnstructuredGrid * const &gridData,
        const std::string &epart_file, const int nElem,
        const FEANode * const &fNode_ptr,
        const IALocal_meshSize * const &lmsize_ptr,
        const IAExtractor * const &ext_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        const BernsteinBasis_Array * const &bern_x,
        const BernsteinBasis_Array * const &bern_y,
        const BernsteinBasis_Array * const &bern_z,
        const double * const * const &pointArrays ) const;


    // --------------------------------------------------------------
    // ! build2doutput: assigns to gridData the sampling points, 
    //                  sampling elements, and data that associated with
    //                  them (similar to the 3d counterpart).
    //                  Note, only 2 BernsteinBasis_Array objects are
    //                  used as input.
    // --------------------------------------------------------------
    void build2doutput(vtkUnstructuredGrid * const &gridData,
        const std::string &epart_file, const int nElem,
        const FEANode * const &fNode_ptr,
        const IALocal_meshSize * const &lmsize_ptr,
        const IAExtractor * const &ext_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        const BernsteinBasis_Array * const &bern_x,
        const BernsteinBasis_Array * const &bern_y,
        const double * const * const &pointArrays ) const;
   

    // --------------------------------------------------------------
    // ! interpolateNURBS: returns the value of the solution at a 
    //                     samping point.
    // \para inputVal: the vector with length nLocBas that gives the 
    //                 solution coefficient within the element;
    // \para elem: the finite element class that includes the 
    //             evaluation of basis functions at sampling points
    // \para output: the double array that returns the evaluation of
    //               the inputVal's solution at sampling points
    // users are responsible for allocating and deleting double arrays.
    // --------------------------------------------------------------
    void interpolateNURBS(const double * const &inputVal,
        const FEAElement * const &elem, double * const &output ) const;
   

    // --------------------------------------------------------------
    // ! function overloading for simultaneously interpolating
    //   two solutions. Use for 2D xy coordinates interpolation 
    // --------------------------------------------------------------
    void interpolateNURBS(const double * const &inputVal_1,
        const double * const &inputVal_2,
        const FEAElement * const &elem, 
        double * const &output_1,
        double * const &output_2 ) const;


    // --------------------------------------------------------------
    // ! function overloading, useful for calculating xyz coordinates
    // --------------------------------------------------------------
    void interpolateNURBS(const double * const &inputVal_1,
        const double * const &inputVal_2,
        const double * const &inputVal_3,
        const FEAElement * const &elem, 
        double * const &output_1,
        double * const &output_2,
        double * const &output_3 ) const;

    
    // --------------------------------------------------------------
    // ! interpolatePts: given the control points of the element, 
    //   evaluate the coordinates of the sampling points and insert
    //   them into vtkPoints
    // \para ctrlPts_x(_y, _z): the nLocBas control points in elem
    // \para elem: the pointer to the element basis functions
    // \vtkpts: the pointer to input point xyz coordinates 
    // --------------------------------------------------------------
    void interpolatePts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts ) const;
   
    
    // --------------------------------------------------------------
    // overloading function for 2D case 
    // --------------------------------------------------------------
    void interpolatePts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts ) const;
    
    
    // --------------------------------------------------------------
    // ! interpolateData: given the inputData, i.e., the solution 
    //   vector we prepared, interpolate it using NURBS, and store
    //   them into vtkData
    //   Note: inputData has the following format: 
    //         [ u11, u21, ..., usize1, u12, u22, ..., usize_n ]
    // \para size: the inputData's dimension, 1 means a simple scalar
    //             type data
    // \para ptoffset: offset
    // \para inputData: the double array for input
    // \para elem: this element,
    // \para vtkData: the array in vtk that we store the interpolated
    //                solution
    // --------------------------------------------------------------
    void interpolateData( const int &size, const int &ptoffset,
        const double * const &inputData,
        const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData ) const;


    // --------------------------------------------------------------
    // ! setHexelem: build the Hexahedron element connectivity with
    //               vtkPoints. In each element, there are
    //               nqp = segs * segt * segu sampling points
    //               these points are indexed as:
    //               ii_s + ii_t * segs + ii_u * seg_s * seg_t.
    // \para segs: number of sampling points in s direction
    // \para segt: number of sampling points in t direction
    // \para segu: number of sampling points in u direction
    // \para ptoffset: offset for this element
    // \gridData: the object that stores the connectivity
    void setHexelem( const int &segs, const int &segt, const int &segu,
        const int &ptoffset, vtkUnstructuredGrid * gridData ) const;


    // --------------------------------------------------------------
    // ! setQuadelem: build Quad element connectivity with vtkPoints
    //                in element, there are nqp = segs * segt elemen
    //                these points are indexed as:
    //                ii_s + ii_t * seg_s
    // \para segs: number of sampling points in s direction
    // \para segt: number of sampling points in t direction
    // \para ptoffset: the offset for this element
    // \para gridData: the object that stores the connectivity
    void setQuadelem( const int &segs, const int &segt,
        const int &ptoffset, vtkUnstructuredGrid * gridData ) const;


    // --------------------------------------------------------------
    // ! read_epart: read the element partition file, usually named as
    //               "epart.h5"
    // \para epart_file: epart.h5, i.e., the element partition from 
    //                   analysis code's preprocessor
    // \para esize: the total element, should be read from global info
    // \para elem_part: the int array with size esize, is the "part" in
    //                  "epart.h5"
    void read_epart( const std::string &epart_file, const int &esize,
       int * const &elem_part ) const;
};

#endif
