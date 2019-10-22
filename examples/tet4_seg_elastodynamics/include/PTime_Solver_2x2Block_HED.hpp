#ifndef PTIME_SOLVER_2X2BLOCK_HED_HPP
#define PTIME_SOLVER_2X2BLOCK_HED_HPP
// ==================================================================
// PTime_Solver_2x2Block_HED.hpp
//
// Parallel time solver for 2x2 Block Hyper-Elasto-Dynamics.
//
// Author: Ju Liu
// Date: Feb. 22 2018
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_Solver_2x2Block_HED.hpp"

class PTime_Solver_2x2Block_HED
{
  public:
    PTime_Solver_2x2Block_HED( 
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const std::string &input_name_d,
        const std::string &input_name_p,
        const std::string &input_name_v,
        const int &input_record_freq, 
        const int &input_renew_tang_freq,
        const double &input_final_time );

    ~PTime_Solver_2x2Block_HED();

    
    void print_info() const;


    void TM_GenAlpha_Solve(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_dot_disp,
        const PDNSolution * const &init_dot_pres,
        const PDNSolution * const &init_dot_velo,
        const PDNSolution * const &init_sol_disp,
        const PDNSolution * const &init_sol_pres,
        const PDNSolution * const &init_sol_velo,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem_2x2Block * const &lassem_ptr,
        IPGAssem_2x2Block * const &gassem_ptr,
        IPLinear_Solver_2x2Block * const &lsolver_ptr,
        PNonlinear_Solver_2x2Block_HED * const &nsolver_ptr );


  private:
    const double final_time;
    const int sol_record_freq;
    const int renew_tang_freq;
    const std::string pb_d_name;
    const std::string pb_p_name;
    const std::string pb_v_name;

    PDNSolution * pre_dot_disp, * pre_dot_pres, * pre_dot_velo;
    PDNSolution * pre_sol_disp, * pre_sol_pres, * pre_sol_velo;
    PDNSolution * cur_dot_disp, * cur_dot_pres, * cur_dot_velo;
    PDNSolution * cur_sol_disp, * cur_sol_pres, * cur_sol_velo;

    std::string Name_Generator( const std::string &base_name,
        const int &counter ) const;
};

#endif
