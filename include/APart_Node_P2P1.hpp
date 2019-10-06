#ifndef APART_NODE_P2P1_HPP
#define APART_NODE_P2P1_HPP
// ==================================================================
// APart_Node_P2P1.hpp
//
// This is the code that generate the local node, the ghost node
// based on the P2 discretization.
// 
// The Gmsh code generate the 10-node tetrahedron with the first 4
// node as the vertex node (meaning the LIEN[ee][0-4] gives the vertex
// node indices.)
//
// We will collect the vertex nodes and create a parallel layout for
// the pressure nodes.
//
// Author: Ju Liu
// Date: Feb. 15 2018
// ==================================================================
#include "APart_Node.hpp"
#include "ALocal_IEN.hpp"

class APart_Node_P2P1 : public APart_Node
{
  public:
    APart_Node_P2P1( const std::string &fileBaseName,
       const int &rank, const ALocal_IEN * const &pien_p2 );

    virtual ~APart_Node_P2P1();

    virtual void print_info() const;
};

#endif
