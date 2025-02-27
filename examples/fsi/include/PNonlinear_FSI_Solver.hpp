#ifndef PNONLINEAR_FSI_SOLVER_HPP
#define PNONLINEAR_FSI_SOLVER_HPP
// ============================================================================
// PNonlinear_FSI_Solver.hpp
//
// This is the nonlienar solver for FSI problems.
//
// Author: Ju Liu
// Date: Jan 4 2022
// ============================================================================
#include "IFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_V.hpp"
#include "PDNSolution_P.hpp"

class PNonlinear_FSI_Solver
{
  public:
    PNonlinear_FSI_Solver(
        std::unique_ptr<IPGAssem> in_gassem_mesh,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_mesh,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<Matrix_PETSc> in_bc_mesh_mat,        
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        std::unique_ptr<IFlowRate> in_flrate,
        std::unique_ptr<PDNSolution> in_sol_base,
        std::unique_ptr<APart_Node> in_pnode_v,
        const double &input_nrtol, const double &input_natol, 
        const double &input_ndtol, const int &input_max_iteration, 
        const int &input_renew_freq, 
        const int &input_renew_threshold );

    PNonlinear_FSI_Solver(
        std::unique_ptr<IPGAssem> in_gassem_prestress,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        std::unique_ptr<APart_Node> in_pnode_v,
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol,
        const int &input_max_iteration, 
        const int &input_renew_freq,
        const int &input_renew_threshold );

    ~PNonlinear_FSI_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    void print_lsolver_mesh_info() const {lsolver_mesh->print_info();}

    void write_prestress_hdf5() const
    {
      SYS_T::print_fatal_if( gassem_prestress == nullptr, "Error: gassem_prestress is nullptr!\n" );
      gassem_prestress->write_prestress_hdf5();
    }

    void GenAlpha_Seg_solve_FSI(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_dot_pres,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_pres,
        const ALocal_InflowBC * const &infnbc_part,
        const IGenBC * const &gbc,
        IPGAssem * const &gassem_ptr,
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &dot_pres,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        PDNSolution * const &pres,
        bool &conv_flag, int &nl_counter ) const;

    void GenAlpha_Seg_solve_Prestress(
        const bool &new_tangent_flag,
        const double &prestress_tol,
        const double &curr_time,
        const double &dt,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_dot_pres,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_pres,
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &dot_pres,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        PDNSolution * const &pres,
        bool &prestress_conv_flag, int &nl_counter ) const;
  
  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    const std::unique_ptr<IPGAssem> gassem_mesh;
    const std::unique_ptr<IPGAssem> gassem_prestress;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver_mesh;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<Matrix_PETSc> bc_mesh_mat;
    const std::unique_ptr<TimeMethod_GenAlpha> tmga;
    const std::unique_ptr<IFlowRate> flrate;
    const std::unique_ptr<PDNSolution> sol_base;
    const std::unique_ptr<const APart_Node> pnode_v;

    void Print_convergence_info( const int &count, const double &rel_err,
        const double &abs_err ) const
    {
      SYS_T::commPrint( "  === NR ite: %d, r_error: %e, a_error: %e \n", count, rel_err, abs_err);
    }

    // This function will add input * val to the solution vector in output, and
    // the ghost values in output will be updated.
    // User is responsible for making sure that the layouts of input and
    // output->solution are identical.
    void update_solid_kinematics( const double &val,
        const Vec &input, PDNSolution * const &output ) const;

    void rescale_inflow_value( const double &stime,
        const ALocal_InflowBC * const &infbc,
        PDNSolution * const &sol ) const;
};

#endif
