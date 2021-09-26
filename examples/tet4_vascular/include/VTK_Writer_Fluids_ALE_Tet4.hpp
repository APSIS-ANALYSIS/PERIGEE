#ifndef VTK_WRITER_FLUIDS_ALE_TET4_HPP
#define VTK_WRITER_FLUIDS_ALE_TET4_HPP
// ==================================================================
// VTK_Writer_Fluids_ALE_Tet4.hpp
// 
// This is a class that specifically designed for the visualization
// for fluid mechanics in ALE coordinates using 4-node tetrahedral 
// element.
//
// This writer is specifically designed for the 7dof formulation.
//
// Date Created: Jan. 25 2017
// ==================================================================
#include "ALocal_Elem.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"
#include "vtkCellData.h"

class VTK_Writer_Fluids_ALE_Tet4
{
  public:
    VTK_Writer_Fluids_ALE_Tet4( const int &in_nelem,
        const std::string &epart_file );

    ~VTK_Writer_Fluids_ALE_Tet4();
    
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
    const int nLocBas;
    const int nElem;

    std::vector<int> epart_map;
};

#endif
