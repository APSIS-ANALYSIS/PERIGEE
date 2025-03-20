#ifndef PTIME_NS_HERK_SOLVER_HPP
#define PTIME_NS_HERK_SOLVER_HPP
// ==================================================================
// PTime_NS_HERK_Solver.hpp
//
// Parallel time solver using HERK for NS.
//
// Author: Yujie Sun
// Date: Mar 10 2025
// ==================================================================
#include "IFlowRate.hpp"
#include "ITimeMethod_RungeKutta.hpp"
#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNTimeStep.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_V.hpp"
#include "PDNSolution_P.hpp"
#include "ALocal_InflowBC.hpp"

class PTime_NS_HERK_Solver
{
  public:
    PTime_NS_HERK_Solver(
        std::unique_ptr<PGAssem_Block_NS_FEM_HERK> in_gassem,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<ITimeMethod_RungeKutta> in_tmRK,
        std::unique_ptr<IFlowRate> in_flrate,
        std::unique_ptr<IFlowRate> in_dot_flrate,
        std::unique_ptr<PDNSolution> in_sol_base,
        std::unique_ptr<ALocal_InflowBC> in_infnbc,
        const std::string &input_name, const int &in_nlocalnode, 
        const int &input_record_freq, const double &input_final_time );

    ~PTime_NS_HERK_Solver() = default;

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    void TM_NS_HERK(
        const bool &restart_init_assembly_flag,
        std::unique_ptr<PDNSolution> init_sol,
        std::unique_ptr<PDNSolution> init_velo,
        std::unique_ptr<PDNSolution> init_dot_velo,
        std::unique_ptr<PDNSolution> init_pres,
        std::unique_ptr<PDNTimeStep> time_info,
        Mat &shell ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const std::string pb_name; // the problem base name for the solution
    const int nlocalnode; // the mumber of nodes (without ghost nodes) in the local area

    const std::unique_ptr<PGAssem_Block_NS_FEM_HERK> gassem;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<ITimeMethod_RungeKutta> tmRK;
    const std::unique_ptr<IFlowRate> flrate;
    const std::unique_ptr<IFlowRate> dot_flrate;
    const std::unique_ptr<PDNSolution> sol_base;
    const std::unique_ptr<const ALocal_InflowBC> infnbc;
  
    std::string Name_Generator( const int &counter ) const;

    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;

    void HERK_Solve_NS(
        const double &curr_time, const double &dt,
        PDNSolution ** const &cur_velo_sols,
        PDNSolution * const &cur_velo,
        PDNSolution * const &cur_dot_velo,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution * const &cur_pres,
        PDNSolution ** const &pre_velo_sols,
        PDNSolution * const &pre_velo,
        PDNSolution ** const &pre_pres_sols,
        PDNSolution * const &pre_pres,
        PDNSolution * const &pre_velo_before,
        PDNSolution * const &cur_sol,
        Mat &shell ) const;

      void rescale_inflow_velo( const double &stime,
          PDNSolution * const &velo ) const;
        
      void rescale_inflow_dot_velo( const double &stime,
          PDNSolution * const &dot_velo ) const;

      void Update_dot_step(const Vec &vp, 
          PDNSolution * const &step) const;

      void Update_pressure_velocity(     
          PDNSolution * const &velo,
          PDNSolution * const &pres,
          const PDNSolution * const &step) const;
  
      void Update_solutions(   
          const PDNSolution * const &velo,
          const PDNSolution * const &pres,
          PDNSolution * const &sol) const;
};

#endif
