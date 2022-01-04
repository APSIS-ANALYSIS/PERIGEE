#include "PNonlinear_FSI_Solver.hpp"

PNonlinear_FSI_Solver::PNonlinear_FSI_Solver(
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol, const int &input_max_iteration,
    const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{}

PNonlinear_FSI_Solver::~PNonlinear_FSI_Solver()
{}

void PNonlinear_FSI_Solver::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::print_sep_line();
}

void PNonlinear_FSI_Solver::rescale_inflow_value( const double &stime,
    const ALocal_Inflow_NodalBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infbc -> get_Num_LD( nbc_id );

    const double val = flrate -> get_flow_rate( nbc_id, stime );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc -> get_LDN( nbc_id, ii );

      const int base_idx[3] = { node_index * 3, node_index * 3 + 1, node_index * 3 + 2 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      VecSetValue(sol->solution, node_index*3  , base_vals[0] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*3+1, base_vals[1] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*3+2, base_vals[2] * val, INSERT_VALUES);
    }
  }

  sol -> Assembly_GhostUpdate();
}

// EOF
