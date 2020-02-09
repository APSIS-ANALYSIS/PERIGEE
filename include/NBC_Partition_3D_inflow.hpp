#ifndef NBC_PARTITION_3D_INFLOW_HPP
#define NBC_PARTITION_3D_INFLOW_HPP
// ==================================================================
// NBC_Partition_3D_inflow.hpp
//
// Inflow Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_inflow, which contains two additional information:
// inflow surface area and the inflow surface outward normal vector.
//
// Date crated: Aug. 9 2017
// Author: Ju Liu
// ==================================================================
#include "NBC_Partition_3D.hpp"

class NBC_Partition_3D_inflow : public NBC_Partition_3D
{
  public:
    NBC_Partition_3D_inflow(const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_3D_inflow();

    virtual void write_hdf5( const char * FileName ) const;

  private:
    double actarea, facearea;

    std::vector<double> outvec;

    int num_out_bc_pts;

    std::vector<double> centroid;

    std::vector<double> outline_pts;
};

#endif
