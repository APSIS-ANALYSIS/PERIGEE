#ifndef PDNSOLUTION_SMOOTH_STRESS_HPP
#define PDNSOLUTION_SMOOTH_STRESS_HPP
// ============================================================================
// PDNSolution_Smooth_Stress.hpp
//
// This is the solution for the solution smoother.
//
// Date: Jan. 17 2024
// ============================================================================
#include "PDNSolution.hpp"

class PDNSolution_Smooth_Stress : public PDNSolution
{
  public:
    PDNSolution_Smooth_Stress( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true,
        const std::string &in_name = "solution_smooth" );

    virtual ~PDNSolution_Smooth_Stress() = default;

  private:
    const std::string sol_name;
    const bool is_print;

    // case 0 : generate a full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );
};

#endif