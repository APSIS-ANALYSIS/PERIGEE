#ifndef PTIME_CMM_SOLVER_HPP
#define PTIME_CMM_SOLVER_HPP
// ==================================================================
// PTime_CMM_Solver.hpp
//
// Parallel time solver for CMM.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_CMM_Solver.hpp"

class PTime_CMM_Solver
{
  public:
    PTime_CMM_Solver( const std::string &input_name, 
        const int &input_record_freq, const int &input_renew_tang_freq, 
        const double &input_final_time );

    ~PTime_CMM_Solver();

    void print_info() const;

    // **** PRESTRESS TODO: additional arg ALocal_Wall_Prestress
    void TM_CMM_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &init_dot_sol,
        const PDNSolution * const &init_sol,
        const PDNSolution * const &init_dot_sol_wall_disp,
        const PDNSolution * const &init_sol_wall_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Inflow_NodalBC * const &infnbc_part,
        const ALocal_Ring_NodalBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
        IGenBC * const &gbc,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_CMM_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter,
        const std::string &prefix = "" ) const;
    
    std::string Name_dot_Generator( const int &counter,
        const std::string &prefix = "" ) const;
    
    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
};

#endif
