#ifndef VTK_WRITER_FSI_HPP
#define VTK_WRITER_FSI_HPP
// ============================================================================
// VTK_Writer_FSI.hpp
//
// This is a class designed for the visualization of FSI problems based on
// linear tet element and trilinear hex element.
//
// Date: Jan 17 2022
// ============================================================================
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"
#include "IMaterialModel.hpp"
#include "Tissue_property.hpp"

#include "vtkIntArray.h"
#include "vtkCellData.h"

class VTK_Writer_FSI
{
  public:
    VTK_Writer_FSI( const int &in_nelem,
        const int &in_nlocbas, 
        const std::string &epart_file );

    ~VTK_Writer_FSI();

    // Write the fluid and solid domain together in a unified continuum body
    void writeOutput(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    void writeOutput_fluid(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const std::vector<int> &fien,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML );

    void writeOutput_solid_cur(
    	const FEANode * const &fnode_ptr,
    	const ALocal_IEN * const &lien_v,
    	const ALocal_IEN * const &lien_p,
    	const std::vector<int> &sien,
    	const ALocal_Elem * const &lelem_ptr,
    	const IVisDataPrep * const &vdata_ptr,
	IMaterialModel ** const &matmodel,
    	FEAElement * const &elemptr,
    	const IQuadPts * const &quad,
	const Tissue_property * const &tp_ptr,
    	const double * const * const &pointArrays,
    	const int &rank, const int &size,
    	const int &num_of_nodes,
    	const double &sol_time,
    	const std::string &outputBName,
    	const std::string &outputName,
    	const bool &isXML );
   
    void writeOutput_solid_ref(
    	const FEANode * const &fnode_ptr,
    	const ALocal_IEN * const &lien_v,
    	const ALocal_IEN * const &lien_p,
    	const std::vector<int> &sien,
    	const ALocal_Elem * const &lelem_ptr,
    	const IVisDataPrep * const &vdata_ptr,
    	FEAElement * const &elemptr,
    	const IQuadPts * const &quad,
    	const double * const * const &pointArrays,
    	const int &rank, const int &size,
    	const int &num_of_nodes,
    	const double &sol_time,
    	const std::string &outputBName,
    	const std::string &outputName,
    	const bool &isXML );

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

    // --------------------------------------------------------------
    // Interpolate von-Mises stress at sampling points
    // --------------------------------------------------------------
    void interpolateVonStress( const int * const &ptid,
        const double * const &ctrlPts_x,
        const double * const &ctrlPts_y,
        const double * const &ctrlPts_z,
	const std::vector<Vector_3> &eleBasis_r,
        const std::vector<Vector_3> &eleBasis_c,
        const std::vector<Vector_3> &eleBasis_l,
        const std::vector<double> &inputDisp,
        const std::vector<double> &inputPres,
        const FEAElement * const &elem,
        IMaterialModel * const &model,
        vtkDoubleArray * const &vtkData );
};

#endif
