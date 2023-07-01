#ifndef VTK_WRITER_FSI_TET4_HPP
#define VTK_WRITER_FSI_TET4_HPP
// ==================================================================
// VTK_Writer_FSI_Tet4.hpp
// 
// This is a class that specifically designed for the visualization
// for FSI problems using 4-node tetrahedral element.
//
// This writer is specifically designed for the VMS-mixed formulation,
// and it will treat FSI problems as a whole unified continuum.
//
// Date Created: Aug. 11 2017
// ==================================================================
#include "ALocal_Elem.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"
#include "vtkCellData.h"

class VTK_Writer_FSI_Tet4
{
  public:
    VTK_Writer_FSI_Tet4( const int &in_nelem,
        const std::string &epart_file );

    ~VTK_Writer_FSI_Tet4();

    // Write the fluid and solid domain together in a unified continuum
    // body    
    void writeOutput(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    // Write only the fluid sub-domain
    void writeOutput_fluid(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const std::vector<int> &fien,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    // Write only the solid sub-domain in current configuration
    void writeOutput_solid(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const std::vector<int> &fien,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    // Write only the solid sub-domain in ref configuration
    void writeOutput_solid_ref(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const std::vector<int> &fien,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

  private:
    const int nLocBas;
    const int nElem;

    std::vector<int> epart_map;

    // --------------------------------------------------------------
    // Interpolate det(F) at sampling points 
    // --------------------------------------------------------------
    void interpolateJ( const int * const &ptid, 
        const std::vector<double> &inputData,
        const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );
};

#endif
