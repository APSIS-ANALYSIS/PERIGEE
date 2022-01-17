#ifndef VTK_WRITER_FSI_TET4_HPP
#define VTK_WRITER_FSI_TET4_HPP
// ============================================================================
// VTK_Writer_FSI_Tet4.hpp
//
// This is a class designed for the visualization of FSI problems based on
// linear tet element.
//
// Date: Jan 17 2022
// ============================================================================
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

  private:
    const int nLocBas, nElem;

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
