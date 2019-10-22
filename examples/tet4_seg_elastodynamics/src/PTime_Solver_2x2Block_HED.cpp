#include "PTime_Solver_2x2Block_HED.hpp"

PTime_Solver_2x2Block_HED::PTime_Solver_2x2Block_HED( 
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const std::string &input_name_d,
    const std::string &input_name_p,
    const std::string &input_name_v,
    const int &input_record_freq,
    const int &input_renew_tang_freq,
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), 
  pb_d_name(input_name_d),
  pb_p_name(input_name_p),
  pb_v_name(input_name_v)
{
  pre_dot_disp = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  pre_dot_pres = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
  pre_dot_velo = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  
  pre_sol_disp = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  pre_sol_pres = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
  pre_sol_velo = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  
  cur_dot_disp = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  cur_dot_pres = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
  cur_dot_velo = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  
  cur_sol_disp = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  cur_sol_pres = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
  cur_sol_velo = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
}



PTime_Solver_2x2Block_HED::~PTime_Solver_2x2Block_HED()
{
  delete pre_dot_disp; pre_dot_disp = NULL;
  delete pre_dot_pres; pre_dot_pres = NULL;
  delete pre_dot_velo; pre_dot_velo = NULL;

  delete pre_sol_disp; pre_sol_disp = NULL;
  delete pre_sol_pres; pre_sol_pres = NULL;
  delete pre_sol_velo; pre_sol_velo = NULL;

  delete cur_dot_disp; cur_dot_disp = NULL;
  delete cur_dot_pres; cur_dot_pres = NULL;
  delete cur_dot_velo; cur_dot_velo = NULL;
  
  delete cur_sol_disp; cur_sol_disp = NULL;
  delete cur_sol_pres; cur_sol_pres = NULL;
  delete cur_sol_velo; cur_sol_velo = NULL;
}


std::string PTime_Solver_2x2Block_HED::Name_Generator( 
    const std::string &base_name, const int &counter ) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name( base_name );
  out_name.append(temp.str());
  return out_name;
}


void PTime_Solver_2x2Block_HED::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "final time: %e \n", final_time);
  PetscPrintf(PETSC_COMM_WORLD, "solution record frequency : %d \n", sol_record_freq);
  PetscPrintf(PETSC_COMM_WORLD, "tangent update frequency over time steps: %d \n", renew_tang_freq);
  PetscPrintf(PETSC_COMM_WORLD, "solution base name: %s %s %s \n", pb_d_name.c_str(),
      pb_p_name.c_str(), pb_v_name.c_str() );
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PTime_Solver_2x2Block_HED::TM_GenAlpha_Solve(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_dot_disp,
    const PDNSolution * const &init_dot_pres,
    const PDNSolution * const &init_dot_velo,
    const PDNSolution * const &init_sol_disp,
    const PDNSolution * const &init_sol_pres,
    const PDNSolution * const &init_sol_velo,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem_2x2Block * const &lassem_ptr,
    IPGAssem_2x2Block * const &gassem_ptr,
    IPLinear_Solver_2x2Block * const &lsolver_ptr,
    PNonlinear_Solver_2x2Block_HED * const &nsolver_ptr )
{
  pre_dot_disp -> Copy( init_dot_disp );
  pre_dot_pres -> Copy( init_dot_pres );
  pre_dot_velo -> Copy( init_dot_velo );

  pre_sol_disp -> Copy( init_sol_disp );
  pre_sol_pres -> Copy( init_sol_pres );
  pre_sol_velo -> Copy( init_sol_velo );

  cur_dot_disp -> Copy( init_dot_disp );
  cur_dot_pres -> Copy( init_dot_pres );
  cur_dot_velo -> Copy( init_dot_velo );

  cur_sol_disp -> Copy( init_sol_disp );
  cur_sol_pres -> Copy( init_sol_pres );
  cur_sol_velo -> Copy( init_sol_velo );

  std::string sol_d_name = Name_Generator( pb_d_name, time_info->get_index() );
  std::string sol_p_name = Name_Generator( pb_p_name, time_info->get_index() );
  std::string sol_v_name = Name_Generator( pb_v_name, time_info->get_index() );

  cur_sol_disp -> WriteBinary( sol_d_name.c_str() );
  cur_sol_pres -> WriteBinary( sol_p_name.c_str() );
  cur_sol_velo -> WriteBinary( sol_v_name.c_str() );

  bool conv_flag, renew_flag;
  int nl_counter;

  bool rest_flag = restart_init_assembly_flag;

  while( time_info -> get_time() < final_time )
  {
    if( time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    nsolver_ptr -> GenAlpha_solve( renew_flag, time_info->get_time(),
        time_info->get_step(), pre_dot_disp, pre_dot_pres, pre_dot_velo, 
        pre_sol_disp, pre_sol_pres, pre_sol_velo, tmga_ptr,
        lien_ptr, feanode_ptr, nbc_part, ebc_part,
        elementv, elements, quad_v, quad_s, lassem_ptr,
        gassem_ptr, lsolver_ptr, cur_dot_disp, cur_dot_pres, cur_dot_velo,
        cur_sol_disp, cur_sol_pres, cur_sol_velo, conv_flag, nl_counter );

    time_info -> TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index());

    if( time_info->get_index()%sol_record_freq == 0)
    {
      sol_d_name = Name_Generator( pb_d_name, time_info->get_index() );
      sol_p_name = Name_Generator( pb_p_name, time_info->get_index() );
      sol_v_name = Name_Generator( pb_v_name, time_info->get_index() );

      cur_sol_disp -> WriteBinary( sol_d_name.c_str() );
      cur_sol_pres -> WriteBinary( sol_p_name.c_str() );
      cur_sol_velo -> WriteBinary( sol_v_name.c_str() );
    }

    pre_dot_disp -> Copy( cur_dot_disp );
    pre_dot_pres -> Copy( cur_dot_pres );
    pre_dot_velo -> Copy( cur_dot_velo );

    pre_sol_disp -> Copy( cur_sol_disp );
    pre_sol_pres -> Copy( cur_sol_pres );
    pre_sol_velo -> Copy( cur_sol_velo );
  }
}

// EOF
