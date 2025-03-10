#ifndef PDNSOLUTION_P_HPP
#define PDNSOLUTION_P_HPP
// ============================================================================
// PDNSolution_P.hpp
//
// This is a solution class for the pressure field, which is attached to a
// pressure-type mesh with 1 dof per grid point.
//
// Author: Ju Liu
// Date: Dec. 29 2021
// ============================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"

class PDNSolution_P : public PDNSolution
{
  public:
    PDNSolution_P( const APart_Node * const &pNode,
        const int &type, const bool &isprint = false,
        const std::string &in_name = "solution_pressure" );

    virtual ~PDNSolution_P() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // --------------------------------------------------------------
    // case 0: generate full zero vector
    // --------------------------------------------------------------
    void Init_zero( const APart_Node * const &pNode );
};

#endif
