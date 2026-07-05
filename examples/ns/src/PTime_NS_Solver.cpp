#include "PTime_NS_Solver.hpp"

void PTime_NS_Solver::write_transport_solution(
    const PDNSolution * const &transport_sol,
    const std::string &file_name,
    const PGAssem_Transport_NS_FEM * const &transport_gassem ) const
{
  if(transport_gassem->has_hi_model())
  {
    PDNSolution hi_sol(*transport_sol);
    HI_T::convert_D_to_HI(&hi_sol, transport_gassem->get_hi_beta());
    hi_sol.WriteBinary(file_name);
  }
  else
  {
    transport_sol->WriteBinary(file_name);
  }
}

PTime_NS_Solver::PTime_NS_Solver(
    std::unique_ptr<PNonlinear_NS_Solver> in_nsolver,
    const std::string &input_name,
    const int &input_record_freq, const int &input_renew_tang_freq,
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name), nsolver(std::move(in_nsolver))
{}

void PTime_NS_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PTime_NS_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
    const std::string &solname ) const
{
  std::ofstream restart_file("restart_file.txt", std::ofstream::out | std::ofstream::trunc);
  if( restart_file.is_open() )
  {
    restart_file<<timeinfo->get_index()<<std::endl;
    restart_file<<timeinfo->get_time()<<std::endl;
    restart_file<<timeinfo->get_step()<<std::endl;
    restart_file<<solname.c_str()<<std::endl;
    restart_file.close();
  }
  else
    SYS_T::print_fatal("Error: PTimeSolver cannot open restart_file.txt");
}

void PTime_NS_Solver::TM_NS_GenAlpha(
    const bool &restart_init_assembly_flag,
    std::unique_ptr<PDNSolution> init_dot_sol,
    std::unique_ptr<PDNSolution> init_sol,
    std::unique_ptr<PDNTimeStep> time_info,
    const ALocal_InflowBC * const &infnbc_part,
    IGenBC * const &gbc,
    IPGAssem * const &gassem_ptr ) const
{
  auto pre_sol     = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto cur_sol     = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto pre_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);
  auto cur_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const auto sol_name = Name_Generator(time_info->get_index());
    cur_sol->WriteBinary(sol_name);

    const auto sol_dot_name = Name_dot_Generator(time_info->get_index());
    cur_dot_sol->WriteBinary(sol_dot_name);
  }

  bool renew_flag;
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
    nl_counter = nsolver->GenAlpha_Solve_NS( renew_flag,
        time_info->get_time(), time_info->get_step(), pre_dot_sol.get(),
        pre_sol.get(), cur_dot_sol.get(), cur_sol.get(), infnbc_part,
        gbc, gassem_ptr );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      const auto sol_name = Name_Generator( time_info->get_index() );
      cur_sol->WriteBinary(sol_name);

      const auto sol_dot_name = Name_dot_Generator(time_info->get_index());
      cur_dot_sol->WriteBinary(sol_dot_name);
    }

    // Calculate the flow rate & averaged pressure on all outlets
    record_outlet_data(cur_sol.get(), cur_dot_sol.get(), time_info.get(), gbc, gassem_ptr, false, true);

    // Calcualte the inlet data
    record_inlet_data(cur_sol.get(), time_info.get(), infnbc_part, gassem_ptr, false, true);

    // Prepare for next time step
    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
  }
}

void PTime_NS_Solver::TM_NS_Transport_GenAlpha(
    const bool &restart_init_assembly_flag,
    std::unique_ptr<PDNSolution> init_dot_sol,
    std::unique_ptr<PDNSolution> init_sol,
    std::unique_ptr<PDNSolution> init_dot_transport,
    std::unique_ptr<PDNSolution> init_transport,
    std::unique_ptr<PDNTimeStep> time_info,
    const ALocal_InflowBC * const &infnbc_part,
    IGenBC * const &gbc,
    IPGAssem * const &gassem_ptr,
    PGAssem_Transport_NS_FEM * const &transport_gassem,
    PLinear_Solver_PETSc * const &transport_lsolver,
    const std::string &transport_name,
    const double &fluid_density ) const
{
  auto pre_sol = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto cur_sol = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto pre_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);
  auto cur_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);

  auto pre_transport = SYS_T::make_unique<PDNSolution>(*init_transport);
  auto cur_transport = SYS_T::make_unique<PDNSolution>(*init_transport);
  auto pre_dot_transport = SYS_T::make_unique<PDNSolution>(*init_dot_transport);
  auto cur_dot_transport = SYS_T::make_unique<PDNSolution>(*init_dot_transport);

  if(restart_init_assembly_flag == false)
  {
    cur_sol->WriteBinary(Name_Generator(time_info->get_index()));
    cur_dot_sol->WriteBinary(Name_dot_Generator(time_info->get_index()));
    write_transport_solution(cur_transport.get(),
        Transport_Name_Generator(transport_name, time_info->get_index()),
        transport_gassem);
    cur_dot_transport->WriteBinary(Transport_Name_dot_Generator(transport_name, time_info->get_index()));
  }

  bool renew_flag;
  int nl_counter = 0;
  bool rest_flag = restart_init_assembly_flag;

  const double alpha_f = nsolver->get_alpha_f();
  const double alpha_m = nsolver->get_alpha_m();
  const double gamma = nsolver->get_gamma();

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    if( nl_counter == 1 ) renew_flag = false;

    nl_counter = nsolver->GenAlpha_Solve_NS( renew_flag,
        time_info->get_time(), time_info->get_step(), pre_dot_sol.get(),
        pre_sol.get(), cur_dot_sol.get(), cur_sol.get(), infnbc_part,
        gbc, gassem_ptr );

    // Predictor for the scalar transport equation.
    cur_transport->Copy(*pre_transport);
    cur_dot_transport->Copy(*pre_dot_transport);
    cur_dot_transport->ScaleValue((gamma - 1.0) / gamma);

    PDNSolution dot_transport_alpha(*pre_dot_transport);
    dot_transport_alpha.ScaleValue(1.0 - alpha_m);
    dot_transport_alpha.PlusAX(*cur_dot_transport, alpha_m);

    PDNSolution transport_alpha(*pre_transport);
    transport_alpha.ScaleValue(1.0 - alpha_f);
    transport_alpha.PlusAX(*cur_transport, alpha_f);

    PDNSolution flow_alpha(*pre_sol);
    flow_alpha.ScaleValue(1.0 - alpha_f);
    flow_alpha.PlusAX(*cur_sol, alpha_f);

    transport_gassem->Clear_KG();
    transport_gassem->Assem_tangent_residual(
        &dot_transport_alpha, &transport_alpha, &flow_alpha,
        time_info->get_time(), time_info->get_step() );

    transport_lsolver->SetOperator(transport_gassem->K);

    PDNSolution transport_step(pre_transport.get());
    transport_lsolver->Solve(transport_gassem->G, &transport_step);

    cur_dot_transport->PlusAX(&transport_step, -1.0);
    cur_transport->PlusAX(&transport_step, -gamma * time_info->get_step());

    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    if( time_info->get_index()%sol_record_freq == 0 )
    {
      cur_sol->WriteBinary(Name_Generator(time_info->get_index()));
      cur_dot_sol->WriteBinary(Name_dot_Generator(time_info->get_index()));
      write_transport_solution(cur_transport.get(),
          Transport_Name_Generator(transport_name, time_info->get_index()),
          transport_gassem);
      cur_dot_transport->WriteBinary(Transport_Name_dot_Generator(transport_name, time_info->get_index()));
    }

    record_outlet_data(cur_sol.get(), cur_dot_sol.get(), time_info.get(), gbc, gassem_ptr, false, true);
    record_inlet_data(cur_sol.get(), time_info.get(), infnbc_part, gassem_ptr, false, true);
    record_outlet_HI_data(cur_sol.get(), cur_transport.get(), time_info.get(), gbc,
        gassem_ptr, transport_gassem, fluid_density, false, true);

    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
    pre_transport->Copy(*cur_transport);
    pre_dot_transport->Copy(*cur_dot_transport);
  }
}

void PTime_NS_Solver::record_outlet_HI_data(
    const PDNSolution * const &flow_sol,
    const PDNSolution * const &transport_sol,
    const PDNTimeStep * const &time_info,
    IGenBC * const &gbc,
    const IPGAssem * const &gassem_ptr,
    const PGAssem_Transport_NS_FEM * const &transport_gassem,
    const double &fluid_density,
    bool is_driver,
    bool is_restart ) const
{
  if(transport_gassem == nullptr || !transport_gassem->has_hi_model()) return;

  auto mode = is_restart ? std::ofstream::app : std::ofstream::trunc;

  for(int ff=0; ff<gbc->get_num_ebc(); ++ff)
  {
    const double flowrate = gassem_ptr->Assem_surface_flowrate(flow_sol, ff);
    const double hi_area_int = transport_gassem->Assem_surface_HI_integral(transport_sol, ff);
    const double hi_flux_int = fluid_density
        * transport_gassem->Assem_surface_HI_flowrate(transport_sol, flow_sol, ff);
    const double hi_ave = std::abs(flowrate) > 1.0e-14
        ? hi_flux_int / (fluid_density * flowrate) : 0.0;

    if( SYS_T::get_MPI_rank() == 0 )
    {
      std::ofstream ofile;
      ofile.open(Outlet_HI_Name_Generator(ff).c_str(), std::ofstream::out | mode);

      if( !is_driver )
      {
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'
             <<flowrate<<'\t'<<hi_area_int<<'\t'<<hi_flux_int<<'\t'<<hi_ave<<std::endl;
      }
      else if( !is_restart )
      {
        ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"Flow rate"<<'\t'
             <<"HI surface integral"<<'\t'<<"HI flux integral"<<'\t'
             <<"Flow-weighted HI"<<'\n';
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'
             <<flowrate<<'\t'<<hi_area_int<<'\t'<<hi_flux_int<<'\t'<<hi_ave<<std::endl;
      }

      ofile.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
}

void PTime_NS_Solver::record_inlet_data(
    const PDNSolution * const &sol,
    const PDNTimeStep * const &time_info,
    const ALocal_InflowBC * const &infnbc_part,
    const IPGAssem * const &gassem_ptr,
    bool is_driver,
    bool is_restart ) const
{
  auto mode = is_restart ? std::ofstream::app : std::ofstream::trunc;

  for(int ff=0; ff<infnbc_part->get_num_nbc(); ++ff)
  {
    const double flrate = gassem_ptr->Assem_surface_flowrate(sol, infnbc_part, ff);

    const double avepre = gassem_ptr->Assem_surface_ave_pressure(sol, infnbc_part, ff);

    if( SYS_T::get_MPI_rank() == 0 )
    {
      std::ofstream ofile;
      ofile.open( gen_flowfile_name("Inlet_", ff).c_str(), std::ofstream::out | mode );

      if( !is_driver )
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'
           <<flrate<<'\t'<<avepre<<std::endl;
      else
      {
        if( !is_restart )
        {
          ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"Flow rate"<<'\t'<<"Face averaged pressure"<<'\n';

          ofile<<time_info->get_index()<<'\t'
             <<time_info->get_time()<<'\t'
             <<flrate<<'\t'<<avepre<<std::endl;
        }
      }

      ofile.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
}

void PTime_NS_Solver::record_outlet_data(
    const PDNSolution * const &sol,
    const PDNSolution * const &dot_sol,
    const PDNTimeStep * const &time_info,
    IGenBC * const &gbc,
    const IPGAssem * const &gassem_ptr,
    bool is_driver,
    bool is_restart) const
{
  auto mode = is_restart ? std::ofstream::app : std::ofstream::trunc;

  for(int ff=0; ff<gbc->get_num_ebc(); ++ff)
  {
    const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate(
        dot_sol, ff);

    const double face_flrate = gassem_ptr -> Assem_surface_flowrate(
        sol, ff);

    const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
        sol, ff);

    double lpn_pressure;

    if ( is_driver )
    {
      gbc -> reset_initial_sol( ff, face_flrate, face_avepre, time_info->get_time(), is_restart );

      lpn_pressure = gbc -> get_P( ff, dot_face_flrate, face_flrate,
        time_info -> get_time() );
    }
    else
    {
      lpn_pressure = gbc -> get_P( ff, dot_face_flrate, face_flrate,
        time_info -> get_time() );

      gbc -> reset_initial_sol( ff, face_flrate, lpn_pressure, time_info->get_time(), false );
    }

    if( SYS_T::get_MPI_rank() == 0 )
    {
      std::ofstream ofile;
      ofile.open( gen_flowfile_name("Outlet_", ff).c_str(), std::ofstream::out | mode );

      if( !is_driver )
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'
             <<dot_face_flrate<<'\t'<<face_flrate<<'\t'
             <<face_avepre<<'\t'<<lpn_pressure<<std::endl;
      else
      {
        if( !is_restart )
        {
          ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"dot Flow rate"<<'\t'<<"Flow rate"<<'\t'
               <<"Face averaged pressure"<<'\t'<<"Reduced model pressure"<<'\n';

          ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'
               <<dot_face_flrate<<'\t'<<face_flrate<<'\t'
               <<face_avepre<<'\t'<<lpn_pressure<<std::endl;
        }
      }

      ofile.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
}

// EOF
