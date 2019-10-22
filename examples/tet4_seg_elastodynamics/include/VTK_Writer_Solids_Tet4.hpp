#ifndef VTK_WRITER_SOLIDS_TET4_HPP
#define VTK_WRITER_SOLIDS_TET4_HPP
// ==================================================================
// VTK_Writer_Solids_Tet4.hpp
// 
// This is a class that specifically designed for the visualization
// for solid mechanics problems using 4-node tetrahedral element.
//
// This writer is specifically designed for the VMS-mixed formulation.
//
// Date Created: Jan. 25 2017
// ==================================================================
#include "ALocal_Elem.hpp"
#include "IMaterialModel.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"
#include "vtkCellData.h"

class VTK_Writer_Solids_Tet4
{
  public:
    VTK_Writer_Solids_Tet4( const int &in_nelem,
        const std::string &epart_file,
        const int &in_nLocBas = 4 );

    ~VTK_Writer_Solids_Tet4();

    
    // Visualize output in the current configuration for isotropic
    // materials: Disp, J, pressure, and velocity 
    void writeOutput_cur_isotropic(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        IMaterialModel * const &model,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    
    // Visualize output in the current configuration for anisotropic
    // materials: Disp, J, average J, Cauchy stress in z direction,
    //            fibre orientation, pressure, and velocity 
    void writeOutput_cur_anisotropic(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        IMaterialModel * const &model,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    
    // Visualize output in the material configuratoin: Disp, Pres, Velo. 
    void writeOutput_ref(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        IMaterialModel * const &model,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );


  private:
    const int nLocBas;
    const int nElem;

    Interpolater intep;

    std::vector<int> epart_map;

    int * IEN_e;
    double * ectrl_x;
    double * ectrl_y;
    double * ectrl_z;

    // --------------------------------------------------------------
    // Interpolate det(F) at sampling points 
    // --------------------------------------------------------------
    void interpolateJ( const int &ptoffset, const double * const &inputData,
        const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );

    // --------------------------------------------------------------
    // Interpolate det(F) and average over element 
    // --------------------------------------------------------------
    void interpolateJele( const int &ptoffset,
        const double * const &inputData, FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );
  
    // --------------------------------------------------------------
    // Interpolate deviatoric Cauchy stress magnitude in the z-direction
    // --------------------------------------------------------------
    void interpolateCauchyZ( const int &ptoffset,
        const double * const &inputData, FEAElement * const &elem,
        IMaterialModel * const &model,
        vtkDoubleArray * const &vtkData );
 
    // --------------------------------------------------------------
    // Interpolate Cauchy stress (with pressure) in the z-direction
    // --------------------------------------------------------------
    void interpolateCauchyZpPres( const int &ptoffset,
        const double * const &inputData_u, 
        const double * const &inputData_p, 
        FEAElement * const &elem,
        IMaterialModel * const &model,
        vtkDoubleArray * const &vtkData );
 
    // --------------------------------------------------------------
    // Interpolate orientation measure := Fa dot Fb
    // --------------------------------------------------------------
    void interpolateOrientation( const int &ptoffset,
        const double * const &inputData, FEAElement * const &elem,
        IMaterialModel * const &model,
        vtkDoubleArray * const &vtkData );
    
    // --------------------------------------------------------------
    // Interpolate the von Mise stress
    // --------------------------------------------------------------
    void interpolateVonMise( const int &ptoffset,
        const double * const &inputData, FEAElement * const &elem,
        IMaterialModel * const &model,
        vtkDoubleArray * const &vtkData );
};


#endif
