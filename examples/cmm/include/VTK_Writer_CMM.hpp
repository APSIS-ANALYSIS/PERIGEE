#ifndef VTK_WRITER_CMM_HPP
#define VTK_WRITER_CMM_HPP
// ==================================================================
// VTK_Writer_CMM.hpp
// 
// This is a class that specifically designed for the visualization
// for fluid mechanics. 
//
// Author: Ju Liu
// Date Created: Feb. 12 2020
// ==================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "ALocal_Elem.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"

class VTK_Writer_CMM
{
  public:
    VTK_Writer_CMM( const int &in_nelem, const int &in_nlocbas, 
        const std::string &epart_file );

    ~VTK_Writer_CMM();
    
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

    Interpolater intep;

    std::vector<int> epart_map;

    std::vector<int> IEN_e;
    
    std::vector<double> ectrl_x, ectrl_y, ectrl_z;
};

#endif
