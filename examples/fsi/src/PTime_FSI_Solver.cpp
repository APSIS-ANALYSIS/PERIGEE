#include "PTime_FSI_Solver.hpp"

PTime_FSI_Solver::PTime_FSI_Solver( const std::string &input_name, 
    const int &input_record_freq, const int &input_renew_tang_freq, 
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
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
    const PDNSolution * const &sol_base,
    const PDNSolution * const &init_dot_disp,
    const PDNSolution * const &init_dot_velo,
    const PDNSolution * const &init_dot_pres,
    const PDNSolution * const &init_disp,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_pres,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_InflowBC * const &infnbc,
    const ALocal_NBC * const &nbc_mesh,
    const ALocal_EBC * const &ebc_v,
    const ALocal_EBC * const &ebc_p,
    const ALocal_EBC * const &ebc_mesh,
    IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const Tissue_prestress * const &ps_ptr,
    IPLocAssem_2x2Block * const &lassem_fluid_ptr,
    IPLocAssem_2x2Block ** const &lassem_solid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
    const PNonlinear_FSI_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * pre_dot_velo = new PDNSolution( init_dot_velo );
  PDNSolution * pre_dot_pres = new PDNSolution( init_dot_pres );

  PDNSolution * pre_disp = new PDNSolution( init_disp );
  PDNSolution * pre_velo = new PDNSolution( init_velo );
  PDNSolution * pre_pres = new PDNSolution( init_pres );

  PDNSolution * cur_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * cur_dot_velo = new PDNSolution( init_dot_velo );
  PDNSolution * cur_dot_pres = new PDNSolution( init_dot_pres );

  PDNSolution * cur_disp = new PDNSolution( init_disp );
  PDNSolution * cur_velo = new PDNSolution( init_velo );
  PDNSolution * cur_pres = new PDNSolution( init_pres );

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
    nsolver_ptr->GenAlpha_Seg_solve_FSI( renew_flag, time_info->get_time(),
        time_info->get_step(), is_v, is_p, sol_base, 
        pre_dot_disp, pre_dot_velo, pre_dot_pres, pre_disp, pre_velo, pre_pres, 
        tmga_ptr, flr_ptr, alelem_ptr, lien_v, lien_p, feanode_ptr, pnode_v, pnode_p, 
        nbc_v, nbc_p, infnbc,
        nbc_mesh, ebc_v, ebc_mesh, gbc, bc_mat, bc_mesh_mat,
        elementv, elements, quad_v, quad_s, ps_ptr, 
        lassem_fluid_ptr, lassem_solid_ptr, lassem_mesh_ptr,
        gassem_ptr, gassem_mesh_ptr, lsolver_ptr, lsolver_mesh_ptr,
        cur_dot_disp, cur_dot_velo, cur_dot_pres, cur_disp, cur_velo, cur_pres, 
        conv_flag, nl_counter );

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
    for(int face=0; face<ebc_v -> get_num_ebc(); ++face)
    {
      // Calculate 3D dot flow rate on the outlets
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_disp,
          cur_dot_velo, lassem_fluid_ptr, elements, quad_s, ebc_v, face );

      // Calculate 3D flow rate on the outlets
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( cur_disp,
          cur_velo, lassem_fluid_ptr, elements, quad_s, ebc_v, face);

      // Calculate 3D averaged pressure on outlets
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp, cur_pres, lassem_fluid_ptr, elements, quad_s, ebc_v, ebc_p, face);

      // Calculate 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate, time_info->get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( ebc_v->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
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
          cur_disp, cur_velo, lassem_fluid_ptr, elements, quad_s, infnbc, face );

      const double inlet_face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_disp, cur_pres, lassem_fluid_ptr, elements, quad_s, infnbc, face );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( infnbc->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    pre_dot_disp -> Copy( cur_dot_disp );
    pre_dot_velo -> Copy( cur_dot_velo );
    pre_dot_pres -> Copy( cur_dot_pres );

    pre_disp -> Copy( cur_disp );
    pre_velo -> Copy( cur_velo );
    pre_pres -> Copy( cur_pres );
  }

  delete pre_dot_disp; delete pre_dot_velo; delete pre_dot_pres;
  delete pre_disp; delete pre_velo; delete pre_pres;

  delete cur_dot_disp; delete cur_dot_velo; delete cur_dot_pres;
  delete cur_disp; delete cur_velo; delete cur_pres;
}

void PTime_FSI_Solver::TM_FSI_Prestress(
    const bool &is_record_sol_flag,
    const double &prestress_tol,
    const IS &is_v,
    const IS &is_p,
    const PDNSolution * const &init_dot_disp,
    const PDNSolution * const &init_dot_velo,
    const PDNSolution * const &init_dot_pres,
    const PDNSolution * const &init_disp,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_pres,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_v,
    const ALocal_NBC * const &nbc_p,
    const ALocal_EBC * const &ebc_v,
    const ALocal_EBC * const &ebc_p,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    Tissue_prestress * const &ps_ptr,
    IPLocAssem_2x2Block ** const &lassem_solid_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    const PNonlinear_FSI_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * pre_dot_velo = new PDNSolution( init_dot_velo );
  PDNSolution * pre_dot_pres = new PDNSolution( init_dot_pres );

  PDNSolution * pre_disp = new PDNSolution( init_disp );
  PDNSolution * pre_velo = new PDNSolution( init_velo );
  PDNSolution * pre_pres = new PDNSolution( init_pres );

  PDNSolution * cur_dot_disp = new PDNSolution( init_dot_disp );
  PDNSolution * cur_dot_velo = new PDNSolution( init_dot_velo );
  PDNSolution * cur_dot_pres = new PDNSolution( init_dot_pres );

  PDNSolution * cur_disp = new PDNSolution( init_disp );
  PDNSolution * cur_velo = new PDNSolution( init_velo );
  PDNSolution * cur_pres = new PDNSolution( init_pres );

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

    // Nullify the solid solutions
    Nullify_solid_dof( pnode_v, 3, pre_dot_disp );
    Nullify_solid_dof( pnode_v, 3, pre_dot_velo );
    Nullify_solid_dof( pnode_p, 1, pre_dot_pres );

    Nullify_solid_dof( pnode_v, 3, pre_disp );
    Nullify_solid_dof( pnode_v, 3, pre_velo );
    Nullify_solid_dof( pnode_p, 1, pre_pres );

    Nullify_solid_dof( pnode_v, 3, cur_dot_disp );
    Nullify_solid_dof( pnode_v, 3, cur_dot_velo );
    Nullify_solid_dof( pnode_p, 1, cur_dot_pres );

    Nullify_solid_dof( pnode_v, 3, cur_disp );
    Nullify_solid_dof( pnode_v, 3, cur_velo );
    Nullify_solid_dof( pnode_p, 1, cur_pres );

    nsolver_ptr -> GenAlpha_Seg_solve_Prestress( renew_flag, prestress_tol,
        time_info->get_time(), time_info->get_step(), is_v, is_p,
        pre_dot_disp, pre_dot_velo, pre_dot_pres, pre_disp, pre_velo, pre_pres,
        tmga_ptr, alelem_ptr, lien_v, lien_p, feanode_ptr, pnode_v, pnode_p,
        nbc_v, nbc_p, ebc_v, ebc_p, bc_mat, elementv, elements, quad_v, quad_s,
        ps_ptr, lassem_solid_ptr, gassem_ptr, lsolver_ptr,
        cur_dot_disp, cur_dot_velo, cur_dot_pres, cur_disp, cur_velo, cur_pres,
        prestress_conv_flag, nl_counter );

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

    pre_dot_disp -> Copy( cur_dot_disp );
    pre_dot_velo -> Copy( cur_dot_velo );
    pre_dot_pres -> Copy( cur_dot_pres );

    pre_disp -> Copy( cur_disp );
    pre_velo -> Copy( cur_velo );
    pre_pres -> Copy( cur_pres );
  }

  delete pre_dot_disp; delete pre_dot_velo; delete pre_dot_pres;
  delete pre_disp; delete pre_velo; delete pre_pres;

  delete cur_dot_disp; delete cur_dot_velo; delete cur_dot_pres;
  delete cur_disp; delete cur_velo; delete cur_pres;
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
