#ifndef VTK_WRITER_SOLID_HPP
#define VTK_WRITER_SOLID_HPP
// ============================================================================
// VTK_Writer_Solid.hpp
//
// Visualization writer for hyperelastic solid problems.
//
// Date: Feb. 01 2026
// ============================================================================
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "FEANode.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "Tensor2_3D.hpp"
#include "SymmTensor2_3D.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"

class VTK_Writer_Solid
{
  public:
    VTK_Writer_Solid( const int &in_nelem,
        const int &in_nlocbas, const std::string &epart_file,
        std::unique_ptr<MaterialModel_Mixed_Elasticity> in_matmodel );

    ~VTK_Writer_Solid() = default;

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
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML,
        const bool &is_ref );

  private:
    const int nLocBas, nElem;
    std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel;

    std::vector<int> epart_map;

    void interpolateJ_Cauchy( const int * const &ptid,
        const double * const &dispData,
        const double * const &presData,
        const FEAElement * const &elem,
        vtkDoubleArray * const &detFData,
        vtkDoubleArray * const &cauchyData ) const;
};

#endif
