#ifndef ALOCAL_IEN_P2P1_HPP
#define ALOCAL_IEN_P2P1_HPP
// ==================================================================
// ALocal_IEN_P2P1.hpp
//
// Local IEN array for the P1 part of the P2/P1 Taylor-Hood element.
//
// This class should be called after the APart_Node_P2P1 has been 
// created. The APart_Node_P2P1 class stores the P1 vertex indices
// using the original P2 layout. In this LIEN array, we will adopt
// the LIEN from the P2 configuration and reset the value based on
// the APart_Node_P2P1 class (since the node list has been re-ordered
// by many operations).
//
// Author: Ju Liu
// Date: Feb. 15 2018
// ==================================================================
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"

class ALocal_IEN_P2P1 : public ALocal_IEN
{
  public:
    ALocal_IEN_P2P1( const std::string &fileBaseName, 
        const int &cpu_rank,
        const APart_Node * const &node_u,
        const APart_Node * const &node_p );

    virtual ~ALocal_IEN_P2P1();

    virtual void print_info() const;

};

#endif
