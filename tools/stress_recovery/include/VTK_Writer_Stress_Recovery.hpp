#ifndef VTK_WRITER_STRESS_RECOVERY_HPP
#define VTK_WRITER_STRESS_RECOVERY_HPP
// ==================================================================
// VTK_Writer_Stress_Recovery.hpp
// 
// This is a class that specifically designed for the visualization
// for smoothed solutions. 
//
// Date Created: Jan. 22 2024
// ==================================================================
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"

class VTK_Writer_Stress_Recovery
{
  public:
    VTK_Writer_Stress_Recovery( const int &in_nelem, const int &in_nlocbas, 
        const std::string &epart_file );

    ~VTK_Writer_Stress_Recovery() = default;
    
    void writeOutput(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

  private:
    const int nLocBas, nElem;

    std::vector<int> epart_map;
};

#endif