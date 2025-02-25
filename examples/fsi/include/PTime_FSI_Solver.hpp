#ifndef PTIME_FSI_SOLVER_HPP
#define PTIME_FSI_SOLVER_HPP
// ============================================================================
// PTime_FSI_Solver.hpp
//
// Parallel time solver for time integration.
//
// Date: Jan. 13 2022
// Author: Ju Liu
// ============================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_FSI_Solver.hpp"

class PTime_FSI_Solver
{
  public:
    PTime_FSI_Solver(
        std::unique_ptr<PNonlinear_NS_Solver> in_nsolver,
        const std::string &input_name,      
        const int &input_record_freq, const int &input_renew_tang_freq, 
        const double &input_final_time );        

    ~PTime_FSI_Solver();

    void print_info() const;

    // ------------------------------------------------------------------------
    // Generate a file name for inlet/outlet face as prefix_xxx_data.txt
    // ------------------------------------------------------------------------
    std::string gen_flowfile_name(const std::string &prefix, const int &id) const
    {
      std::ostringstream ss;
      ss << prefix;

      if(id < 10) ss << "00";
      else if(id < 100) ss << "0";

      ss << id << "_data.txt";
      return ss.str();
    }

    void TM_FSI_GenAlpha(
        const bool &restart_init_assembly_flag,
        const IS &is_v,
        const IS &is_p,
        std::unique_ptr<PDNSolution> init_dot_disp,
        std::unique_ptr<PDNSolution> init_dot_velo,
        std::unique_ptr<PDNSolution> init_dot_pres,
        std::unique_ptr<PDNSolution> init_disp,
        std::unique_ptr<PDNSolution> init_velo,
        std::unique_ptr<PDNSolution> init_pres,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const ALocal_InflowBC * const &infnbc,
        IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr,
        IPGAssem * const &gassem_ptr,
        IPGAssem * const &gassem_mesh_ptr ) const;

    void TM_FSI_Prestress(
        const bool &is_record_sol_flag,
        const double &prestress_tol,
        const IS &is_v,
        const IS &is_p,
        std::unique_ptr<PDNSolution> init_dot_disp,
        std::unique_ptr<PDNSolution> init_dot_velo,
        std::unique_ptr<PDNSolution> init_dot_pres,
        std::unique_ptr<PDNSolution> init_disp,
        std::unique_ptr<PDNSolution> init_velo,
        std::unique_ptr<PDNSolution> init_pres,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        Tissue_prestress * const &ps_ptr,
        IPGAssem * const &gassem_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    const std::unique_ptr<PNonlinear_NS_Solver> nsolver;

    std::string Name_Generator( const std::string &middle_name, 
        const int &counter ) const;

    std::string Name_dot_Generator( const std::string &middle_name, 
        const int &counter ) const;

    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;

    void Nullify_solid_dof( const APart_Node * const &pnode,
        const int &in_dof, PDNSolution * const &sol ) const;
};

#endif
