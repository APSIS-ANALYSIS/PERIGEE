#ifndef PDNSOLUTION_ELASTODYNAMICS_HPP
#define PDNSOLUTION_ELASTODYNAMICS_HPP
// ============================================================================
// PDNSolution_Elastodynamics.hpp
//
// This is the solution for elastodynamics problems with 3 dof per node.
//
// Date: Oct 24 2023
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_Elastodynamics : public PDNSolution
{
  public:
    PDNSolution_Elastodynamics( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true,
        const std::string &in_name = "solution_elastodynamics" );

    virtual ~PDNSolution_Elastodynamics();

  private:
    const std::string sol_name;
    const bool is_print;

    // case 0 : generate a full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

};

#endif
