#ifndef PDNSOLUTION_TRANSPORT_HPP
#define PDNSOLUTION_TRANSPORT_HPP
// ============================================================================
// PDNSolution_Transport.hpp
//
// This is the solution for transport problems with 1 dof per node.
//
// Date: Jan 21 2022
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_Transport : public PDNSolution
{
  public:
    PDNSolution_Transport( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true,
        const std::string &in_name = "solution_transport" );

    virtual ~PDNSolution_Transport() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // case 0 : generate a full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );
};

#endif
