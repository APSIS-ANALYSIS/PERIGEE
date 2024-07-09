#ifndef PTIME_NS_SOLVER_HPP
#define PTIME_NS_SOLVER_HPP
// ==================================================================
// PTime_NS_Solver.hpp
//
// Parallel time solver for NS.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_NS_Solver.hpp"

class PTime_NS_Solver
{
  public:
    PTime_NS_Solver( const std::string &input_name, 
        const int &input_record_freq, const int &input_renew_tang_freq, 
        const double &input_final_time );

    ~PTime_NS_Solver() = default;

    void print_info() const;

    void TM_NS_GenAlpha(
        const bool &restart_init_assembly_flag,
        const double &TI_std_dev,
        PDNSolution * const &sol_base,
        const PDNSolution * const &init_dot_sol,
        const PDNSolution * const &init_sol,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ICVFlowRate * const flr_ptr,
        const APart_Node * const &pNode_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_InflowBC * const &infnbc_part,
        const ALocal_EBC * const &ebc_part,
        IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part,
        ALocal_Interface * const &itf_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_NS_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;
    
    std::string Name_dot_Generator( const int &counter ) const;

    std::string Name_disp_Generator( const int &counter ) const;
    
    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
};

#endif
