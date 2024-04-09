#ifndef VTK_WRITER_SMOOTH_FSI_HPP
#define VTK_WRITER_SMOOTH_FSI_HPP
// ==================================================================
// VTK_Writer_Smooth_FSI.hpp
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
#include "IMaterialModel.hpp"
#include "Tissue_property.hpp"
#include "vtkIntArray.h"

class VTK_Writer_Smooth_FSI
{
  public:
    VTK_Writer_Smooth_FSI( const int &in_nelem, const int &in_nlocbas, 
        const std::string &epart_file );

    ~VTK_Writer_Smooth_FSI();
    
    void writeOutput(
      const FEANode * const &fnode_ptr,
    	const ALocal_IEN * const &lien_v,
    	const ALocal_IEN * const &lien_p,
    	const std::vector<int> &sien,
    	const ALocal_Elem * const &lelem_ptr,
    	const IVisDataPrep * const &vdata_ptr,
      IMaterialModel * const &matmodel,
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

  private:
    const int nLocBas, nElem;

    std::vector<int> epart_map;

/*
    void interpolateVonMises( const int * const &ptid,
        const std::vector<double> &inputCauchy,
        const std::vector<double> &inputPres,
        const FEAElement * const &elem,
        vtkDoubleArray * const &vtkData );
*/
    void interpolateF( const int * const &ptid,
      const std::vector<double> &inputGradDisp,
      const FEAElement * const &elem,
      vtkDoubleArray * const &vtkData );
    
    void interpolateStrain( const int * const &ptid,
      const std::vector<double> &inputGradDisp,
      const FEAElement * const &elem,
      IMaterialModel * const &model,
      vtkDoubleArray * const &vtkData );

    void interpolateCauchy( const int * const &ptid,
      const std::vector<Vector_3> &eleBasis_r,
      const std::vector<Vector_3> &eleBasis_c,
      const std::vector<Vector_3> &eleBasis_l,
      const std::vector<double> &inputGradDisp,
      const FEAElement * const &elem,
      IMaterialModel * const &model,
      vtkDoubleArray * const &vtkData );   

    void interpolateVonMises( const int * const &ptid,
      const std::vector<Vector_3> &eleBasis_r,
      const std::vector<Vector_3> &eleBasis_c,
      const std::vector<Vector_3> &eleBasis_l,
      const std::vector<double> &inputGradDisp,
      const std::vector<double> &inputPres,
      const FEAElement * const &elem,
      IMaterialModel * const &model,
      vtkDoubleArray * const &vtkData );    

    void interpolateVonMises_nop( const int * const &ptid,
      const std::vector<Vector_3> &eleBasis_r,
      const std::vector<Vector_3> &eleBasis_c,
      const std::vector<Vector_3> &eleBasis_l,
      const std::vector<double> &inputGradDisp,
      const FEAElement * const &elem,
      IMaterialModel * const &model,
      vtkDoubleArray * const &vtkData );  
};

#endif