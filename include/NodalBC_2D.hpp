#ifndef NODALBC_2D_HPP
#define NODALBC_2D_HPP
// ==================================================================
// NodalBC_2D.hpp
//
// This is an instantiation of INodalBC for 2D problems.
//
// Date: Nov. 17 2017
// Author: Ju Liu
// ==================================================================
#include "INodalBC.hpp"

class NodalBC_2D : public INodalBC
{
  public:
    // Default constructor:
    // Clear the dir_nodes, per_slave_nodes, per_master_nodes;
    // Set ID based on the above "no-nodal bc" setting.
    NodalBC_2D( const int &nFunc );


    // Set dir_nodes as in_pt[ ed_idx ]
    NodalBC_2D( const int &ed_idx,
       const std::vector<std::vector<int> > &in_pt, const int &nFunc );


    // Set dir_nodes as in_pt[ edge_idx[ii] ], for 0<=ii<edge_idx.size().
    NodalBC_2D( const std::vector<int> &edge_idx,
       const std::vector<std::vector<int> > &in_pt, const int &nFunc );


    virtual ~NodalBC_2D();

  private:
    NodalBC_2D() {}; // dis-allow default constructor
};

#endif
