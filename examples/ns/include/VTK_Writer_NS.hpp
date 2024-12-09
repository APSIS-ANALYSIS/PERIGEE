#ifndef VTK_WRITER_NS_HPP
#define VTK_WRITER_NS_HPP
// ==================================================================
// VTK_Writer_NS.hpp
// 
// This is a class that specifically designed for the visualization
// for fluid mechanics. 
//
// Author: Ju Liu
// Date Created: Feb. 12 2020
// ==================================================================
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "FEANode.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"

class VTK_Writer_NS
{
  public:
    VTK_Writer_NS( const int &in_nelem, const int &in_nlocbas, 
        const std::string &epart_file );

    ~VTK_Writer_NS() = default;
    
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
