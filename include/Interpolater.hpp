#ifndef INTERPOLATER_HPP
#define INTERPOLATER_HPP
// ============================================================================
// Interpolater.hpp
// This is a class that designed for repeated evaluation of finite 
// element interpolation. For example, it can be called in the 
// evaluation of sampling points in the visualization routine.
// This class contains a container for R, R_x, R_y, R_z.
//
// To call the interpolate functions, the elements should have already
// been built at sampling points.
//
// Author: Ju Liu
// Date: July 12 2016
// ============================================================================
#include "FEAElement.hpp"
#include "vtkPoints.h"
#include "vtkDoubleArray.h"

class Interpolater
{
  public:
    // ------------------------------------------------------------------------
    // Construct the FE basis containers based on the input
    // nLocBas : number of basis functions and
    // isDer   : the flag telling if the derivatives are needed
    // ------------------------------------------------------------------------
    Interpolater( const int &in_nlocbas );
    
    virtual ~Interpolater();
    
    // print the basic info on screen.
    void print_info() const;

    // ------------------------------------------------------------------------
    // interpolateFE : return a value based on the inputVal, which is
    // the coefficients for the element elem's basis functions.
    // \para inputVal: the vector with length nLocBas that gives the 
    //                 solution coefficient within the element;
    // \para elem: the finite element class that includes the 
    //             evaluation of basis functions at sampling points
    // \para output: the double array that returns the evaluation of
    //               the inputVal's solution at sampling points
    // This function is overloaded to efficiently calculate 1, 2, and
    // 3 fields. This is used, for example, in evaluate physical coor
    // dinates in the visualization routines.
    // ------------------------------------------------------------------------
    void interpolateFE( const double * const &inputVal,
        const FEAElement * const &elem, std::vector<double> &output );

    void interpolateFE( const double * const &inputVal_1,
        const double * const &inputVal_2,
        const FEAElement * const &elem,
        std::vector<double> &output_1, std::vector<double> &output_2 );

    void interpolateFE( const double * const &inputVal_1,
        const double * const &inputVal_2,
        const double * const &inputVal_3,
        const FEAElement * const &elem,
        std::vector<double> &output_1, std::vector<double> &output_2,
        std::vector<double> &output_3 );

    // ------------------------------------------------------------------------
    // interpolateFE_Grad : return a value based on the interpolation
    // of the elements' basis functions' 1st-order derivatives.
    // Output gives the interpolate of dR_dx, dR_dy, dR_dz, respectively
    // ------------------------------------------------------------------------
    void interpolateFE_Grad( const double * const &inputVal,
        const FEAElement * const &elem, std::vector<double> &output_dx,
        std::vector<double> &output_dy, std::vector<double> &output_dz );
    
    void interpolateFE_Grad( const std::vector<double> &inputVal,
        const FEAElement * const &elem, std::vector<double> &output_dx,
        std::vector<double> &output_dy, std::vector<double> &output_dz )
    {interpolateFE_Grad(&inputVal[0], elem, output_dx, output_dy, output_dz);}
    
    // ------------------------------------------------------------------------
    // interpolateVTKPts : given the control points of the element, evaluate
    // the coordinates of the sampling points and insert them into 
    // vtkPoints object which is the output of this funciton.
    // \para ptoffset : offset for the inserting points in vtkPoints object 
    // \para ctrlPts_x/y/z : nLocBas control points
    // \para elem : element pointer
    // \output vtkpts : the pointer to the sampling pts' coordinates
    //
    // Function is overloaded for 2D coordinates, differentiating by the
    // control points' input.
    // ------------------------------------------------------------------------
    void interpolateVTKPts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts );

    void interpolateVTKPts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts );
  
    // ------------------------------------------------------------------------
    // interpolateVTKPts : overloaded for the case of four points with
    //                     given point indices.
    //                     The points id are stored in ptid;
    //                     The number of points are stored in elem by
    //                     elem -> get_numQuapts().
    //                     Users are responsible to make sure ptid.size()
    //                     equals elem -> get_numQuapts().
    // ------------------------------------------------------------------------
    void interpolateVTKPts( const int * const &ptid,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts );

    // ------------------------------------------------------------------------
    // interpolateVTKPts : overloaded for Lagrangian mesh points.
    // The coordinates are updated based on the input displacements.
    // The control points is x = X + u.
    // X : Reference coordinates;
    // u : displacement results;
    // /para disp_vect : the displacement given in the following format 
    // ux_1 uy_1 uz_1 ux_2 uy_2 ... ux_nlocbas ny_nlocbas nz_nlocbas
    // ------------------------------------------------------------------------
    void interpolateVTKPts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const double * const &disp_vect,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts );

    void interpolateVTKPts( const int &ptoffset,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const std::vector<double> &disp_vect,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts )
    {interpolateVTKPts(ptoffset, ctrlPts_x, ctrlPts_y, ctrlPts_z, 
        &disp_vect[0], elem, vtkpts);}

    // ------------------------------------------------------------------------
    // interpolateVTKPts : overloaded for Lagrangian mesh points.
    //                     The points id are stored in ptid;
    //                     The number of points are stored in elem by
    //                     elem -> get_numQuapts().
    //                     Users are responsible to make sure ptid.size()
    //                     equals elem -> get_numQuapts().
    // ------------------------------------------------------------------------
    void interpolateVTKPts( const int * const &ptid,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const double * const &disp_vect,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts );

    void interpolateVTKPts( const int * const &ptid,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
        const std::vector<double> &disp_vect,
        const FEAElement * const &elem,
        vtkPoints * const &vtkpts )
    {interpolateVTKPts(ptid, ctrlPts_x, ctrlPts_y, ctrlPts_z, &disp_vect[0], elem, vtkpts);}

    // ------------------------------------------------------------------------
    // interpolateVTKData : given the inputData, i.e., the solution 
    // vector for the element. Interpolate the solution vector and write
    // them into vtkDoubleArray object as output.
    // \para size: input arrays dimension. If 1, a scalar data is interpolated.
    //             The input data should have the following format:
    //       [ u1_1, u2_1, ..., usize_1, u1_2, u2_2, ..., usize_n ]
    // \para ptoffset
    // \para inputData : the double array for input; or in std::vector
    //                   format, with length nLocBas x size.
    // \para elem : the element pointer
    // \output vtkData: the array in vtk that store the output with 
    //                  the number of component = size
    // ------------------------------------------------------------------------
    void interpolateVTKData( const int &size, const int &ptoffset,
        const double * const &inputData, const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );

    void interpolateVTKData( const int &size, const int &ptoffset,
        const std::vector<double> &inputData, const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData )
    {interpolateVTKData(size, ptoffset, &inputData[0], elem, vtkData);}

    // ------------------------------------------------------------------------
    // interpolateVTKData : overload the previous function by replacing
    //                      automatic numbering by specific given point
    //                      indeices.
    //                      The points id are stored in ptid;
    //                      The number of points are stored in elem by
    //                      elem -> get_numQuapts().
    //                      Users are responsible to make sure ptid.size()
    //                      equals elem -> get_numQuapts().
    // ------------------------------------------------------------------------
    void interpolateVTKData( const int &size, const int * const &ptid, 
        const double * const &inputData, const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );

    void interpolateVTKData( const int &size, const int * const &ptid, 
        const std::vector<double> &inputData, const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData )
    {interpolateVTKData(size, ptid, &inputData[0], elem, vtkData);}

    // ------------------------------------------------------------------------
    // interpolateData : given the input Data, i.e., the element sol-
    //                   ution vector, interpolate the actual value
    //                   and write into a std::vector.
    // \para size: input arrays dimension. If 1, a scalar data is interpolated.
    //             The input data should have the following format:
    //       [ u1_1, u2_1, ..., usize_1, u1_2, u2_2, ..., usize_n ]
    // \para inputData : the double array for input
    // \para elem : the element pointer
    // \output outData: the std::vector that stores the interpolated
    //                  values. Its size is size x nqp. 
    // ------------------------------------------------------------------------
    void interpolateData( const int &size, 
        const double * const &inputData,
        const FEAElement * const &elem, 
        std::vector< std::vector<double> > &outData );

  private:
    // number of local(elemental) basis functions
    const int nLocBas;

    // Disable empty constructor
    Interpolater() = delete;
};

#endif
