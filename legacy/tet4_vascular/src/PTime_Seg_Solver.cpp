#include "PTime_Seg_Solver.hpp"

PTime_Seg_Solver::PTime_Seg_Solver(
    const std::string &input_name, const int &input_record_freq,
    const int &input_renew_tang_freq, const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{}


PTime_Seg_Solver::~PTime_Seg_Solver()
{}


std::string PTime_Seg_Solver::Name_Generator(const int &counter) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}


std::string PTime_Seg_Solver::Name_dot_Generator(const int &counter) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name("dot_");
  out_name.append(pb_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_Seg_Solver::print_info() const
{
  SYS_T::commPrint( "----------------------------------------------------------- \n");
  SYS_T::commPrint( "final time: %e \n", final_time);
  SYS_T::commPrint( "solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint( "tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint( "solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint( "----------------------------------------------------------- \n");
}


void PTime_Seg_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
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


void PTime_Seg_Solver::TM_ALE_NS_GenAlpha( 
    const bool &restart_init_assembly_flag,
    const bool &is_ale_flag,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_NBC * const &nbc_mesh_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_mesh_part,
    IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_fluid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
    PNonlinear_Seg_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  std::string sol_name ("");
  std::string sol_dot_name ("");

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    sol_name = Name_Generator(time_info->get_index());
    cur_disp->WriteBinary(sol_name.c_str());
    
    sol_dot_name = Name_dot_Generator(time_info->get_index());
    cur_velo->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag, renew_flag;
  int nl_counter;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  // Enter into time integration
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    // Inoke the nonlinear equation solver
    nsolver_ptr->GenAlpha_Seg_solve_ALE_NS( renew_flag, is_ale_flag, 
        time_info->get_time(), time_info->get_step(), 
        sol_base, pre_velo, pre_disp, tmga_ptr, flr_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, infnbc_part,
        nbc_mesh_part, ebc_part, ebc_mesh_part, gbc, bc_mat, bc_mesh_mat, 
        elementv, elements, quad_v, quad_s, lassem_fluid_ptr, lassem_mesh_ptr,
        gassem_ptr, gassem_mesh_ptr, lsolver_ptr, lsolver_mesh_ptr,
        cur_velo, cur_disp, conv_flag, nl_counter );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());

      sol_dot_name = Name_dot_Generator(time_info->get_index());
      cur_velo->WriteBinary(sol_dot_name.c_str());
    }

    // Calculate the flow rate & averaged pressure on all outlets
    for(int face=0; face<ebc_part -> get_num_ebc(); ++face)
    {
      // Calculate the 3D dot flow rate on the outlet
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_velo, lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate the 3D flow rate on the outlet
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_disp, lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate the 3D averaged pressure on the outlet
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure( 
          cur_disp, lassem_fluid_ptr, elements, quad_s, ebc_part, face);

      // Calculate the 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate, time_info->get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      // On the CPU 0, write the time, flow rate, averaged pressure, and 0D
      // calculated pressure into the txt file, which is first generated in the
      // driver
      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( ebc_part->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp; delete pre_velo; delete cur_velo;
}


void PTime_Seg_Solver::TM_FSI_GenAlpha( 
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_NBC * const &nbc_mesh_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_mesh_part,
    IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const Tissue_prestress * const &ps_ptr,
    IPLocAssem * const &lassem_fluid_ptr,
    IPLocAssem * const &lassem_solid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
    PNonlinear_Seg_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // Do not overwrite solution if this is a restart 
  if( restart_init_assembly_flag == false )
  {
    const std::string sol_name = Name_Generator(time_info->get_index());
    cur_disp->WriteBinary(sol_name.c_str());

    const std::string sol_dot_name = Name_dot_Generator(time_info->get_index());
    cur_velo->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag, renew_flag;
  int nl_counter;

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

    nsolver_ptr->GenAlpha_Seg_solve_FSI( renew_flag, time_info->get_time(),
        time_info->get_step(), sol_base, pre_velo, pre_disp, tmga_ptr, flr_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, infnbc_part,
        nbc_mesh_part, ebc_part, ebc_mesh_part, gbc, bc_mat, bc_mesh_mat, 
        elementv, elements, quad_v, quad_s, ps_ptr, lassem_fluid_ptr,
        lassem_solid_ptr, lassem_mesh_ptr,
        gassem_ptr, gassem_mesh_ptr, lsolver_ptr, lsolver_mesh_ptr,
        cur_velo, cur_disp, conv_flag, nl_counter );

    time_info->TimeIncrement();

    SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str() );

    if( time_info->get_index()%sol_record_freq == 0)
    {
      const std::string sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());

      const std::string sol_dot_name = Name_dot_Generator(time_info->get_index());
      cur_velo->WriteBinary(sol_dot_name.c_str());
    }

    // Calculate the flow rate on all outlets
    for(int face=0; face<ebc_part -> get_num_ebc(); ++face)
    {
      // Calculate 3D dot flow rate on the outlets
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_velo,
          lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate 3D flow rate on the outlets
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_disp,
          lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate 3D averaged pressure on outlets
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp, lassem_fluid_ptr, elements, quad_s, ebc_part, face);

      // Calculate 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate, time_info->get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( ebc_part->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
        ofile.close();
      }

      MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Write all 0D solutions into a file
    if( SYS_T::get_MPI_rank() == 0 )
      gbc -> write_0D_sol ( time_info->get_index(), time_info->get_time() );

    MPI_Barrier(PETSC_COMM_WORLD);

    // Calculate the flow rate and averaged pressure on all inlets
    for(int face=0; face<infnbc_part -> get_num_nbc(); ++face)
    {
      const double inlet_face_flrate = gassem_ptr -> Assem_surface_flowrate(
          cur_disp, lassem_fluid_ptr, elements, quad_s, infnbc_part, face );

      const double inlet_face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp, lassem_fluid_ptr, elements, quad_s, infnbc_part, face );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( infnbc_part->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Prepare for the next time step
    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp; delete pre_velo; delete cur_velo;
}


void PTime_Seg_Solver::TM_FSI_Prestress(
    const bool &is_record_sol_flag,
    const double &prestress_tol, 
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
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
    Tissue_prestress * const &ps_ptr,
    IPLocAssem * const &lassem_solid_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Seg_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  bool prestress_conv_flag = false, renew_flag;
  int nl_counter = nsolver_ptr -> get_non_max_its();

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

    // nullify the solid solution
    SEG_SOL_T::Insert_zero_solid_UV( anode_ptr, pre_velo ); 
    SEG_SOL_T::Insert_zero_solid_UV( anode_ptr, pre_disp ); 
    SEG_SOL_T::Insert_zero_solid_UV( anode_ptr, cur_velo ); 
    SEG_SOL_T::Insert_zero_solid_UV( anode_ptr, cur_disp ); 
    
    nsolver_ptr->GenAlpha_Solve_Prestress( renew_flag, prestress_tol, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part,
        ebc_part, bc_mat, elementv, elements, quad_v, quad_s, ps_ptr,
        lassem_solid_ptr, gassem_ptr, lsolver_ptr,
        cur_velo, cur_disp, prestress_conv_flag, nl_counter );
    
    time_info->TimeIncrement();

    SYS_T::commPrint( "Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str() );

    if( is_record_sol_flag && time_info->get_index()%sol_record_freq == 0)
    {
      const std::string sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());

      const std::string sol_dot_name = Name_dot_Generator(time_info->get_index());
      cur_velo->WriteBinary(sol_dot_name.c_str());
    }

    // Prepare for the next time step
    pre_disp->Copy(*cur_disp); pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp; delete pre_velo; delete cur_velo;
}

// EOF