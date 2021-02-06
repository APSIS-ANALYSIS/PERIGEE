#ifndef PDNSOLUTION_WALL_DISP_HPP
#define PDNSOLUTION_WALL_DISP_HPP
// ==================================================================
// PDNSolution_Wall_Disp.hpp
//
// This is the solution vector for CMM solver -- wall displacement.
//
// The solution has 3 dofs, that is the three componenet of displacement.
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"

class PDNSolution_Wall_Disp : public PDNSolution
{
  public:
    PDNSolution_Wall_Disp( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const int &type, const bool &isprint = true );

    ~PDNSolution_Wall_Disp();

    // case 0 : generate a full zero solution vector
    void Init_zero( const APart_Node * const &pNode_ptr );

  private:
    const bool is_print;
};

#endif
