#ifndef PDNSOLUTION_V_HPP
#define PDNSOLUTION_V_HPP
// ============================================================================
// PDNSolution_V.hpp
//
// This is a solution class for the kinematics field (displacement, velocity,
// mesh displacement, mesh velocity, etc.), which is attached to a velocity-type
// mesh with 3 dofs per grid point.
//
// Author: Ju Liu
// Date: Dec. 29 2021 
// ============================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_InflowBC.hpp"

class PDNSolution_V : public PDNSolution
{
  public:
    PDNSolution_V( const APart_Node * const &pNode,
        const int &type, const bool &isprint = false,
        const std::string &in_name = "solution_velocity" );

    virtual ~PDNSolution_V() {};

  private:
    const std::string sol_name;
    const bool is_print;

    // --------------------------------------------------------------
    // case 0: generate full zero vector 
    // --------------------------------------------------------------
    void Init_zero( const APart_Node * const &pNode );
};

#endif
