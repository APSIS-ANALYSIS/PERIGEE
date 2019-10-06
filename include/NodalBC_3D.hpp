#ifndef NODALBC_3D_HPP
#define NODALBC_3D_HPP
// ==================================================================
// NodalBC_3D.hpp
//
// This is an instantiation of INodalBC for 3D problems.
//
// Date: Dec. 6 2017
// Author: Ju Liu
// ==================================================================
#include "INodalBC.hpp"

class NodalBC_3D : public INodalBC
{
  public:
    // Constructor: clear the dir_nodes, per_slave_nodes, per_master_nodes, 
    // set ID array based on the "no-nodal-bc" setting.
    NodalBC_3D( const int &nFunc );

    // Set dir_nodes at in_pt[ face_idx ]
    NodalBC_3D( const int &face_idx,
        const std::vector<std::vector<int> > &in_pt, const int &nFunc );

    // Set dir_nodes as in_pt[ face_idx[ii] ], for 0<=ii<face_idx.size().
    NodalBC_3D( const std::vector<int> &face_idx,
        const std::vector<std::vector<int> > &in_pt, const int &nFunc );


    // General NodalBC constructor. The specific implementation of the
    // nodal BC is given from the private functions and are specified
    // from the type flag.
    NodalBC_3D( const std::vector<int> &face_idx,
        const std::vector<std::vector<int> > &in_pt, const int &nFunc,
        const int &type );

    virtual ~NodalBC_3D();

  private:
    NodalBC_3D() {}; // dis-allow default constructor

    // --------------------------------------------------------------    
    // BC_type_1: There are two surfaces given. The first is a Dirichlet
    // surface; the second one has a master-slave mapping, one node
    // is the master while the rest nodes are all slaves.
    // This is used for the tensile test of arterial strips where
    // the traction face is not allowed to deform.
    // --------------------------------------------------------------    
    void BC_type_1( const std::vector<int> &face_idx,
        const std::vector<std::vector<int> > &in_pt, const int &nFunc );
    
    // --------------------------------------------------------------    
    // BC_type_2: We set a master-slave mapping for all faces given
    // in the face_idx list. One node in the face surface is set as 
    // the master node and the rest node are all slaves. This, for
    // example, can be used for tensile test with traction force applied
    // at both ends with face deformation fixed at both ends.
    // --------------------------------------------------------------    
    void BC_type_2( const std::vector<int> &face_idx,
        const std::vector<std::vector<int> > &in_pt, const int &nFunc );
};

#endif
