#include "PTime_NS_Solver.hpp"

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
    std::unique_ptr<PDNSolution> init_mdisp,
    std::unique_ptr<PDNSolution> init_mvelo,
    std::unique_ptr<PDNTimeStep> time_info,
    std::unique_ptr<ALocal_InflowBC> infnbc_part,
    std::unique_ptr<ALocal_RotatedBC> rotnbc_part,
    std::unique_ptr<IGenBC> gbc,
    std::unique_ptr<SI_rotation_info> rot_info,
    std::unique_ptr<IPGAssem> gassem_ptr,
    Mat &shell ) const
{
  PDNSolution * pre_sol = new PDNSolution(*init_sol);
  PDNSolution * cur_sol = new PDNSolution(*init_sol);
  PDNSolution * pre_dot_sol = new PDNSolution(*init_dot_sol);
  PDNSolution * cur_dot_sol = new PDNSolution(*init_dot_sol);
  PDNSolution * pre_disp_mesh = new PDNSolution(*init_mdisp);
  PDNSolution * cur_disp_mesh = new PDNSolution(*init_mdisp);
  PDNSolution * alpha_disp_mesh = new PDNSolution(*init_mdisp);
  PDNSolution * pre_velo_mesh = new PDNSolution(*init_mvelo);
  PDNSolution * cur_velo_mesh = new PDNSolution(*init_mvelo);
  PDNSolution * alpha_velo_mesh = new PDNSolution(*init_mvelo);  

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const auto sol_name = Name_Generator(time_info->get_index());
    cur_sol->WriteBinary(sol_name);
    
    const auto sol_dot_name = Name_dot_Generator(time_info->get_index());
    cur_dot_sol->WriteBinary(sol_dot_name);
  
    const auto sol_disp_name = Name_disp_Generator(time_info->get_index());
    cur_disp_mesh->WriteBinary(sol_disp_name);

    const auto sol_mvelo_name = Name_mvelo_Generator(time_info->get_index());
    cur_velo_mesh->WriteBinary(sol_mvelo_name);
  }

  bool conv_flag, renew_flag;
  int nl_counter = 0;

  bool rest_flag = restart_init_assembly_flag;

  const double alpha_f = nsolver->get_alpha_f();

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  const APart_Node * pNode_ptr = gassem_ptr->Get_pnode();
  const FEANode * feanode_ptr = gassem_ptr->Get_fnode();

  // Enter into time integration
  while( time_info->get_time() < final_time )
  {
    //Calculate the cur_velo_mesh and cur_disp_mesh
    Vec lvelo_mesh, ldisp_mesh, lvelo_alp_mesh, ldisp_alp_mesh;
    double * array_cur_velo_mesh, * array_cur_disp_mesh;
    double * array_alp_velo_mesh, * array_alp_disp_mesh;
    VecGhostGetLocalForm(cur_velo_mesh->solution, &lvelo_mesh);
    VecGhostGetLocalForm(cur_disp_mesh->solution, &ldisp_mesh);
    VecGhostGetLocalForm(alpha_velo_mesh->solution, &lvelo_alp_mesh);
    VecGhostGetLocalForm(alpha_disp_mesh->solution, &ldisp_alp_mesh);    
    VecGetArray(lvelo_mesh, &array_cur_velo_mesh);
    VecGetArray(ldisp_mesh, &array_cur_disp_mesh);
    VecGetArray(lvelo_alp_mesh, &array_alp_velo_mesh);
    VecGetArray(ldisp_alp_mesh, &array_alp_disp_mesh);

    for( int ii=0; ii<pNode_ptr->get_nlocalnode_rotated(); ++ii )
    { 
      // Update the coordinates of the rotated nodes
      const Vector_3 init_pt_xyz = feanode_ptr->get_ctrlPts_xyz(pNode_ptr->get_node_loc_rotated(ii));
      const Vector_3 curr_pt_xyz = get_currPts(init_pt_xyz, time_info->get_time() + time_info->get_step(), rot_info.get()); //get_currPts() may be writtern into Sl_tools
      const Vector_3 aplha_pt_xyz = get_currPts(init_pt_xyz, time_info->get_time() + alpha_f * time_info->get_step(), rot_info.get());

      const Vector_3 radius_alpha = get_radius(aplha_pt_xyz, rot_info.get()); 
      const Vector_3 velo_mesh_alpha = Vec3::cross_product(rot_info->get_angular_velo(time_info->get_time() + alpha_f * time_info->get_step())*rot_info->get_direction_rotated(), radius_alpha);

      const Vector_3 radius_curr = get_radius(curr_pt_xyz, rot_info.get()); //get_radius() may be writtern into Sl_tools  
      const Vector_3 velo_mesh_curr = Vec3::cross_product(rot_info->get_angular_velo(time_info->get_time() + time_info->get_step())*rot_info->get_direction_rotated(), radius_curr);

      const int offset = pNode_ptr->get_node_loc_rotated(ii) * 3;   

      for(int jj=0; jj<3; ++jj)
      {
        array_cur_velo_mesh[offset + jj] = velo_mesh_curr(jj);
        array_cur_disp_mesh[offset + jj] = curr_pt_xyz(jj)-init_pt_xyz(jj);  

        array_alp_velo_mesh[offset + jj] = velo_mesh_alpha(jj);
        array_alp_disp_mesh[offset + jj] = aplha_pt_xyz(jj)-init_pt_xyz(jj);  
      }
    }

    VecRestoreArray(lvelo_mesh, &array_cur_velo_mesh);
    VecRestoreArray(ldisp_mesh, &array_cur_disp_mesh);
    VecRestoreArray(lvelo_alp_mesh, &array_alp_velo_mesh);
    VecRestoreArray(ldisp_alp_mesh, &array_alp_disp_mesh);    
    VecGhostRestoreLocalForm(cur_velo_mesh->solution, &lvelo_mesh);
    VecGhostRestoreLocalForm(cur_disp_mesh->solution, &ldisp_mesh);
    VecGhostRestoreLocalForm(alpha_velo_mesh->solution, &lvelo_alp_mesh);
    VecGhostRestoreLocalForm(alpha_disp_mesh->solution, &ldisp_alp_mesh);

    cur_velo_mesh->GhostUpdate();
    cur_disp_mesh->GhostUpdate();
    alpha_velo_mesh->GhostUpdate();
    alpha_disp_mesh->GhostUpdate();

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
    nsolver->GenAlpha_Solve_NS( renew_flag, 
        time_info->get_time(), time_info->get_step(), 
        pre_dot_sol, pre_sol, pre_velo_mesh, pre_disp_mesh,
        infnbc_part.get(), rotnbc_part.get(), gbc.get(), gassem_ptr.get(),
        cur_dot_sol, cur_sol, cur_velo_mesh, cur_disp_mesh,
        alpha_velo_mesh, alpha_disp_mesh,
        conv_flag, nl_counter, shell );

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

      const auto sol_disp_name = Name_disp_Generator(time_info->get_index());
      cur_disp_mesh->WriteBinary(sol_disp_name);

      const auto sol_mvelo_name = Name_mvelo_Generator(time_info->get_index());
      cur_velo_mesh->WriteBinary(sol_mvelo_name);
    }

    // Calculate the flow rate & averaged pressure on all outlets
    for(int face=0; face<gbc -> get_num_ebc(); ++face)
    {
      // Calculate the 3D dot flow rate on the outlet
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_dot_sol, face); 

      // Calculate the 3D flow rate on the outlet
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_sol, face); 

      // Calculate the 3D averaged pressure on the outlet
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure( 
          cur_sol, face);

      // Calculate the 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate,
         time_info -> get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      // On the CPU 0, write the time, flow rate, averaged pressure, and 0D
      // calculated pressure into the txt file, which is first generated in the
      // driver
      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( gen_flowfile_name("Outlet_", face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<dot_face_flrate<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
   
    // Calcualte the inlet data
    for(int face=0; face<infnbc_part -> get_num_nbc(); ++face)
    {
      const double inlet_face_flrate = gassem_ptr -> Assem_surface_flowrate(
          cur_sol, infnbc_part.get(), face ); 

      const double inlet_face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_sol, infnbc_part.get(), face );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( gen_flowfile_name("Inlet_", face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
        ofile.close();
      } 
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Prepare for next time step
    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
    pre_disp_mesh->Copy(*cur_disp_mesh);    
    pre_velo_mesh->Copy(*cur_velo_mesh);
  }

  delete pre_sol; delete cur_sol; delete pre_dot_sol; delete cur_dot_sol; 
  delete pre_velo_mesh; delete cur_velo_mesh; delete pre_disp_mesh; 
  delete cur_disp_mesh; delete alpha_velo_mesh; delete alpha_disp_mesh;
}

// EOF
