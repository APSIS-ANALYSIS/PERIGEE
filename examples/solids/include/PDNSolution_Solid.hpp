#ifndef PDNSOLUTION_SOLID_HPP
#define PDNSOLUTION_SOLID_HPP
// ============================================================================
// PDNSolution_Solid.hpp
//
// This is a solution class for solid problems with configurable dofs per node.
//
// Date: Feb. 01 2026
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_Solid : public PDNSolution
{
  public:
    PDNSolution_Solid( const APart_Node * const &pNode,
        const int &in_dof_num, const int &type,
        const bool &isprint = false,
        const std::string &in_name = "solution_solid" );

    virtual ~PDNSolution_Solid() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // --------------------------------------------------------------
    // case 0: generate full zero vector
    // --------------------------------------------------------------
    void Init_zero( const APart_Node * const &pNode, const int &in_dof_num );
};

#endif
