#ifndef PDNSOLUTION_V_HPP
#define PDNSOLUTION_V_HPP
// ============================================================================
// PDNSolution_V.hpp
//
// This is a solution class for the velocity field which is attached to a velocity-type
// mesh with 3 dofs per grid point.
//
// Author: Yujie Sun
// Date: Apr. 14 2025 
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_V final : public PDNSolution
{
  public:
    PDNSolution_V( const APart_Node * const &pNode,
        const int &type, const bool &isprint = false,
        const std::string &in_name = "solution_velocity" );

    virtual ~PDNSolution_V() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // --------------------------------------------------------------
    // case 0: generate full zero vector 
    // --------------------------------------------------------------
    void Init_zero( const APart_Node * const &pNode );
};

#endif
