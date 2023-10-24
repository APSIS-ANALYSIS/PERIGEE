#include "PTime_LinearPDE_Solver.hpp"

PTime_LinearPDE_Solver::PTime_LinearPDE_Solver( const std::string &input_name,
        const int &input_record_freq, const int &input_renew_tang_freq,
        const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{}

PTime_LinearPDE_Solver::~PTime_LinearPDE_Solver()
{}

void PTime_LinearPDE_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

std::string PTime_LinearPDE_Solver::Name_Generator( const std::string &middle_name,
    const int &counter ) const
{
  const int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(middle_name);
  out_name.append(temp.str());
  return out_name;
}

std::string PTime_LinearPDE_Solver::Name_dot_Generator( const std::string &middle_name,
    const int &counter ) const
{
  const int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name("dot_");
  out_name.append(pb_name);
  out_name.append(middle_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_LinearPDE_Solver::TM_GenAlpha_Transport(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_dot_sol,
    const PDNSolution * const &init_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_LinearPDE_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_sol = new PDNSolution(*init_sol);
  PDNSolution * cur_sol = new PDNSolution(*init_sol);
  PDNSolution * pre_dot_sol = new PDNSolution(*init_dot_sol);
  PDNSolution * cur_dot_sol = new PDNSolution(*init_dot_sol);

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const std::string sol_name = Name_Generator("temp_", time_info->get_index());
    cur_sol->WriteBinary(sol_name.c_str());

    const std::string sol_dot_name = Name_dot_Generator("temp_", time_info->get_index());
    cur_dot_sol->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag, renew_flag;
  int nl_counter = 0;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());
 
 // Enter into time integration
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    // If the previous step is solved in ONE Newton iteration, we do not update
    // the tangent matrix
    if( nl_counter == 1 ) renew_flag = false;

    // Call the nonlinear equation solver
    nsolver_ptr->GenAlpha_Solve_Transport( renew_flag,
        time_info->get_time(), time_info->get_step(),
        pre_dot_sol, pre_sol, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part,
        ebc_part, bc_mat, elementv, elements, quad_v, quad_s, lassem_ptr,
        gassem_ptr, lsolver_ptr, cur_dot_sol, cur_sol, conv_flag, nl_counter );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      const std::string sol_name = Name_Generator("temp_", time_info->get_index() );
      cur_sol->WriteBinary(sol_name.c_str());

      const std::string sol_dot_name = Name_dot_Generator("temp_", time_info->get_index());
      cur_dot_sol->WriteBinary(sol_dot_name.c_str());
    }

    // Prepare for next time step
    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
  } 

  delete pre_sol; delete cur_sol; delete pre_dot_sol; delete cur_dot_sol;
}

void PTime_LinearPDE_Solver::TM_GenAlpha_Elastodynamics(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_dot_disp,
    const PDNSolution * const &init_dot_velo,
    const PDNSolution * const &init_disp,
    const PDNSolution * const &init_velo,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_LinearPDE_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution( init_disp );
  PDNSolution * cur_disp = new PDNSolution( init_disp );
  PDNSolution * pre_velo = new PDNSolution( init_velo );
  PDNSolution * cur_velo = new PDNSolution( init_velo );
  PDNSolution * pre_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * cur_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * pre_dot_velo = new PDNSolution( init_dot_velo );
  PDNSolution * cur_dot_velo = new PDNSolution( init_dot_velo );

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    std::string sol_name = Name_Generator("disp_", time_info->get_index());
    cur_disp->WriteBinary(sol_name.c_str());

    sol_name = Name_Generator("velo_", time_info->get_index());
    cur_velo->WriteBinary(sol_name.c_str());

    std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
    cur_dot_disp->WriteBinary(sol_dot_name.c_str());

    sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
    cur_dot_velo->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag, renew_flag;
  int nl_counter = 0;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());
 
 // Enter into time integration
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    // If the previous step is solved in ONE Newton iteration, we do not update
    // the tangent matrix
    if( nl_counter == 1 ) renew_flag = false;

    // Call the nonlinear equation solver
    nsolver_ptr->GenAlpha_Solve_Elastodynamics( renew_flag,
        time_info->get_time(), time_info->get_step(),
        pre_dot_disp, pre_dot_velo, pre_disp, pre_velo, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part,
        ebc_part, bc_mat, elementv, elements, quad_v, quad_s, lassem_ptr,
        gassem_ptr, lsolver_ptr, cur_dot_disp, cur_dot_velo,
        cur_disp, cur_velo, conv_flag, nl_counter );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      std::string sol_name = Name_Generator("disp_", time_info->get_index());
      cur_disp->WriteBinary(sol_name.c_str());

      sol_name = Name_Generator("velo_", time_info->get_index());
      cur_velo->WriteBinary(sol_name.c_str());

      std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
      cur_dot_disp->WriteBinary(sol_dot_name.c_str());

      sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
      cur_dot_velo->WriteBinary(sol_dot_name.c_str());
    }

    // Prepare for next time step
    pre_disp -> Copy( cur_disp );
    pre_velo -> Copy( cur_velo );
    pre_dot_disp -> Copy( cur_dot_disp );
    pre_dot_velo -> Copy( cur_dot_velo );
  } 

  delete pre_disp; delete cur_disp; delete pre_dot_disp; delete cur_dot_disp;
  delete pre_velo; delete cur_velo; delete pre_dot_velo; delete cur_dot_velo;
}

// EOF
