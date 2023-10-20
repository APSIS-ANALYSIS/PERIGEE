#ifndef PDNSOLUTION_LinearPDE_HPP
#define PDNSOLUTION_LinearPDE_HPP
// ============================================================================
// PDNSolution_LinearPDE.hpp
//
// This is the solution for linear pde problems with prescribed dof per node.
//
// Date: Jan 21 2022
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_LinearPDE : public PDNSolution
{
  public:
    PDNSolution_LinearPDE( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true,
        const std::string &in_name );

    virtual ~PDNSolution_LinearPDE();

  private:
    const std::string sol_name;
    const bool is_print;

    // case 0 : generate a full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

};

#endif
