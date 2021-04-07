#ifndef NODALBC_3D_RING_HPP
#define NODALBC_3D_RING_HPP
// ==================================================================
// NodalBC_3D_ring.hpp
// 
// This is an instantiation of INodalBC for 3D ring boundary conditions
// enabling in-plane motion for the inlet and outlet ring/boundary
// nodes in CMM.
//
// For ring boundary conditions, 1. there is no periodic type
// boundary condition; 2. For a given inlet/outlet, the ring nodes
// should only be constrained in the dominant component of the
// unit normal.
//
// Author: Ju Liu, Ingrid Lan
// Date: Apr. 7 2021
// ==================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"

class NodalBC_3D_ring : public INodalBC
{
  public:
    // --------------------------------------------------------------
    // Generate an empty ring type boundary condition class
    // --------------------------------------------------------------
    NodalBC_3D_ring(const int &nFunc);
    
    // --------------------------------------------------------------
    // Generate the ring bc given by the inflow and outflow files, 
    // their unit normals, and the wall file.
    // --------------------------------------------------------------
    NodalBC_3D_ring( const std::string &inflow_file,
        const std::vector<double> &inflow_outward_vec,
        const std::string &wallfile,
        const std::vector<std::string> &outflow_files,
        const std::vector< std::vector<double> > &outflow_outward_vec,
        const int &nFunc,
        const int &elemtype = 501 );

    virtual ~NodalBC_3D_ring() {};

  private:
    NodalBC_3D_ring() {};
    
};

#endif
