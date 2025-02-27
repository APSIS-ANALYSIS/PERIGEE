#include "PTime_FSI_Solver.hpp"

PTime_FSI_Solver::PTime_FSI_Solver(
    std::unique_ptr<PNonlinear_FSI_Solver> in_nsolver,
    std::unique_ptr<APart_Node> in_pnode_v,
    std::unique_ptr<APart_Node> in_pnode_p,
    const std::string &input_name,      
    const int &input_record_freq, const int &input_renew_tang_freq, 
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name), 
  nsolver(std::move(in_nsolver)), pnode_v(std::move(in_pnode_v)),
  pnode_p(std::move(in_pnode_p))
{}

PTime_FSI_Solver::~PTime_FSI_Solver()
{}

void PTime_FSI_Solver::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint( "final time: %e \n", final_time);
  SYS_T::commPrint( "solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint( "tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint( "solution base name: %s \n", pb_name.c_str());
  SYS_T::print_sep_line();
}

std::string PTime_FSI_Solver::Name_Generator( const std::string &middle_name,
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

std::string PTime_FSI_Solver::Name_dot_Generator( const std::string &middle_name,
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

void PTime_FSI_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
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

void PTime_FSI_Solver::TM_FSI_GenAlpha(
    const bool &restart_init_assembly_flag,
    const IS &is_v,
    const IS &is_p,
    std::unique_ptr<PDNSolution> init_dot_disp,
    std::unique_ptr<PDNSolution> init_dot_velo,
    std::unique_ptr<PDNSolution> init_dot_pres,
    std::unique_ptr<PDNSolution> init_disp,
    std::unique_ptr<PDNSolution> init_velo,
    std::unique_ptr<PDNSolution> init_pres,
    std::unique_ptr<PDNTimeStep> time_info,
    const ALocal_InflowBC * const &infnbc,
    IGenBC * const &gbc,
    IPGAssem * const &gassem_ptr ) const
{
  auto pre_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto pre_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto pre_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto pre_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto pre_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto pre_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  auto cur_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto cur_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto cur_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto cur_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto cur_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto cur_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  // Do NOT overwrite solution if this is a restart
  if( restart_init_assembly_flag == false )
  {
    std::string sol_name = Name_Generator("disp_", time_info->get_index());
    cur_disp->WriteBinary(sol_name);

    sol_name = Name_Generator("velo_", time_info->get_index());
    cur_velo->WriteBinary(sol_name);

    sol_name = Name_Generator("pres_", time_info->get_index());
    cur_pres->WriteBinary(sol_name);

    std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
    cur_dot_disp->WriteBinary(sol_dot_name);

    sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
    cur_dot_velo->WriteBinary(sol_dot_name);

    sol_dot_name = Name_dot_Generator("pres_", time_info->get_index());
    cur_dot_pres->WriteBinary(sol_dot_name);
  }

  int nl_counter = 100;
  bool renew_flag;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str() );

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    // If the previous step is solved in ONE Newton iteration, we do not update
    // the tangent matrix
    if( nl_counter == 1 ) renew_flag = false;

    bool conv_flag;
    nsolver -> GenAlpha_Seg_solve_FSI( renew_flag, time_info->get_time(),time_info->get_step(), 
        is_v, is_p, pre_dot_disp.get(), pre_dot_velo.get(), pre_dot_pres.get(), pre_disp.get(), 
        pre_velo.get(), pre_pres.get(), infnbc, gbc, gassem_ptr, cur_dot_disp.get(), cur_dot_velo.get(), 
        cur_dot_pres.get(), cur_disp.get(), cur_velo.get(), cur_pres.get(), conv_flag, nl_counter );

    time_info->TimeIncrement();

    SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str() );

    if( time_info->get_index()%sol_record_freq == 0)
    {
      std::string sol_name = Name_Generator("disp_", time_info->get_index());
      cur_disp->WriteBinary(sol_name);

      sol_name = Name_Generator("velo_", time_info->get_index());
      cur_velo->WriteBinary(sol_name);

      sol_name = Name_Generator("pres_", time_info->get_index());
      cur_pres->WriteBinary(sol_name);

      std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
      cur_dot_disp->WriteBinary(sol_dot_name);

      sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
      cur_dot_velo->WriteBinary(sol_dot_name);

      sol_dot_name = Name_dot_Generator("pres_", time_info->get_index());
      cur_dot_pres->WriteBinary(sol_dot_name);
    }

    // Calculate the flow rate on all outlets
    for(int face=0; face<gbc -> get_num_ebc(); ++face)
    {
      // Calculate 3D dot flow rate on the outlets
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_disp.get(),
          cur_dot_velo.get(), face );

      // Calculate 3D flow rate on the outlets
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_disp.get(),
          cur_velo.get(), face);

      // Calculate 3D averaged pressure on outlets
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp.get(), cur_pres.get(), face);

      // Calculate 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate, time_info->get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( gen_flowfile_name("Outlet_", face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
        ofile.close();
      }

      MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Write all 0D solutions into a file
    if( SYS_T::get_MPI_rank() == 0 )
      gbc -> write_0D_sol ( time_info->get_index(), time_info->get_time() );

    // Calculate the flow rate and averaged pressure on all inlets
    for(int face=0; face<infnbc -> get_num_nbc(); ++face)
    {
      const double inlet_face_flrate = gassem_ptr -> Assem_surface_flowrate(
          cur_disp.get(), cur_velo.get(), infnbc, face );

      const double inlet_face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp.get(), cur_pres.get(), infnbc, face );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( gen_flowfile_name("Inlet_", face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    pre_dot_disp -> Copy( cur_dot_disp.get() );
    pre_dot_velo -> Copy( cur_dot_velo.get() );
    pre_dot_pres -> Copy( cur_dot_pres.get() );

    pre_disp -> Copy( cur_disp.get() );
    pre_velo -> Copy( cur_velo.get() );
    pre_pres -> Copy( cur_pres.get() );
  }
}

void PTime_FSI_Solver::TM_FSI_Prestress(
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
    std::unique_ptr<PDNTimeStep> time_info ) const
{
  auto pre_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto pre_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto pre_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto pre_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto pre_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto pre_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  auto cur_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto cur_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto cur_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto cur_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto cur_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto cur_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  bool prestress_conv_flag = false, renew_flag;
  int nl_counter = nsolver -> get_non_max_its();

  SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str() );

  while( time_info->get_time() < final_time && !prestress_conv_flag )
  {
    if(time_info->get_index() % renew_tang_freq == 0) renew_flag = true;
    else renew_flag = false;

    // If the previous step is solved in ONE Newton iteration, we do not update
    // the tangent matrix
    if( nl_counter == 1 ) renew_flag = false;

    // Nullify the solid solutions
    Nullify_solid_dof( pnode_v.get(), 3, pre_dot_disp.get() );
    Nullify_solid_dof( pnode_v.get(), 3, pre_dot_velo.get() );
    Nullify_solid_dof( pnode_p.get(), 1, pre_dot_pres.get() );

    Nullify_solid_dof( pnode_v.get(), 3, pre_disp.get() );
    Nullify_solid_dof( pnode_v.get(), 3, pre_velo.get() );
    Nullify_solid_dof( pnode_p.get(), 1, pre_pres.get() );

    Nullify_solid_dof( pnode_v.get(), 3, cur_dot_disp.get() );
    Nullify_solid_dof( pnode_v.get(), 3, cur_dot_velo.get() );
    Nullify_solid_dof( pnode_p.get(), 1, cur_dot_pres.get() );

    Nullify_solid_dof( pnode_v.get(), 3, cur_disp.get() );
    Nullify_solid_dof( pnode_v.get(), 3, cur_velo.get() );
    Nullify_solid_dof( pnode_p.get(), 1, cur_pres.get() );

    nsolver -> GenAlpha_Seg_solve_Prestress( renew_flag, prestress_tol,
        time_info->get_time(), time_info->get_step(), is_v, is_p, pre_dot_disp.get(), 
        pre_dot_velo.get(), pre_dot_pres.get(), pre_disp.get(), pre_velo.get(), 
        pre_pres.get(), cur_dot_disp.get(), cur_dot_velo.get(), cur_dot_pres.get(), 
        cur_disp.get(), cur_velo.get(), cur_pres.get(), prestress_conv_flag, nl_counter );

    time_info->TimeIncrement();

    SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str() );

    // Record the solution under prestress generation, if meets criteria
    if( is_record_sol_flag && time_info->get_index()%sol_record_freq == 0)
    {
      std::string sol_name = Name_Generator("disp_", time_info->get_index());
      cur_disp->WriteBinary(sol_name);

      sol_name = Name_Generator("velo_", time_info->get_index());
      cur_velo->WriteBinary(sol_name);

      sol_name = Name_Generator("pres_", time_info->get_index());
      cur_pres->WriteBinary(sol_name);

      std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
      cur_dot_disp->WriteBinary(sol_dot_name);

      sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
      cur_dot_velo->WriteBinary(sol_dot_name);

      sol_dot_name = Name_dot_Generator("pres_", time_info->get_index());
      cur_dot_pres->WriteBinary(sol_dot_name);
    }

    pre_dot_disp -> Copy( cur_dot_disp.get() );
    pre_dot_velo -> Copy( cur_dot_velo.get() );
    pre_dot_pres -> Copy( cur_dot_pres.get() );

    pre_disp -> Copy( cur_disp.get() );
    pre_velo -> Copy( cur_velo.get() );
    pre_pres -> Copy( cur_pres.get() );
  }
}

void PTime_FSI_Solver::Nullify_solid_dof( const APart_Node * const &pnode,
    const int &in_dof, PDNSolution * const &sol ) const
{
  SYS_T::print_fatal_if(sol->get_dof_num() != in_dof,
      "Error: PTime_FSI_Solver::Nullify_solid_dof, the input dof value is wrong. \n");

  SYS_T::print_fatal_if(pnode->get_dof() != in_dof,
      "Error: PTime_FSI_Solver::Nullify_solid_dof, the input dof value is wrong. \n");

  SYS_T::print_fatal_if(sol->get_nlocal() != pnode->get_nlocalnode() * in_dof,
      "Error: PTime_FSI_Solver::Nullify_solid_dof, the input solution dimension is wrong. \n");

  Vec lsol;
  VecGhostGetLocalForm(sol->solution, &lsol);

  double * array_sol;
  VecGetArray(lsol, &array_sol);

  const int nlocal_solid = pnode->get_nlocalnode_solid();

  for(int ii=0; ii<nlocal_solid; ++ii)
  {
    const int offset = pnode -> get_node_loc_solid(ii) * in_dof;
    for(int jj=0; jj<in_dof; ++jj)
      array_sol[offset + jj] = 0.0;
  }

  VecRestoreArray(lsol, &array_sol);
  VecGhostRestoreLocalForm(sol->solution, &lsol);

  sol->GhostUpdate(); // update the ghost slots
}

// EOF
