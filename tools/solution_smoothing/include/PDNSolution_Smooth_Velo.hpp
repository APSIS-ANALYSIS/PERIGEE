#ifndef PDNSOLUTION_SMOOTH_Velo_HPP
#define PDNSOLUTION_SMOOTH_Velo_HPP
// ============================================================================
// PDNSolution_Smooth_Vol.hpp
//
// This is the solution for the solution smoother.
//
// Date: Jan. 17 2024
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_Smooth_Velo : public PDNSolution
{
  public:
    PDNSolution_Smooth_Velo( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true,
        const std::string &in_name = "solution_smooth_vol" );

    virtual ~PDNSolution_Smooth_Velo() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // case 0 : generate a full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );
};

#endif