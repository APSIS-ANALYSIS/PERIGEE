#include "PTime_Tet4_NS_3D_Solver.hpp"

PTime_Tet4_NS_3D_Solver::PTime_Tet4_NS_3D_Solver( 
    const std::string &input_name,
    const int &input_record_freq, const int &input_renew_tang_freq,
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{
}

PTime_Tet4_NS_3D_Solver::~PTime_Tet4_NS_3D_Solver()
{}

void PTime_Tet4_NS_3D_Solver::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "final time: %e \n", final_time);
  PetscPrintf(PETSC_COMM_WORLD, "solution record frequency : %d \n", sol_record_freq);
  PetscPrintf(PETSC_COMM_WORLD, "tangent update frequency over time steps: %d \n", renew_tang_freq);
  PetscPrintf(PETSC_COMM_WORLD, "solution base name: %s \n", pb_name.c_str());
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}

std::string PTime_Tet4_NS_3D_Solver::Name_Generator( const int &counter ) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}


void PTime_Tet4_NS_3D_Solver::TM_VMS_GenAlpha(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Tet4_NS_3D_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  std::string sol_name = Name_Generator(time_info->get_index());
  cur_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  int nl_counter;
  bool rest_flag = restart_init_assembly_flag;

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    nsolver_ptr->GenAlpha_VMS_solve( renew_flag, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, ebc_part,
        bc_mat, elementv, elements, quad_v, quad_s, lassem_ptr,
        gassem_ptr, lsolver_ptr, cur_velo, cur_disp, conv_flag, nl_counter );

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index());

    if( time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());
    }

    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }
  delete pre_disp; delete cur_disp; delete pre_velo; delete cur_velo;
}


void PTime_Tet4_NS_3D_Solver::TM_VMS_GenAlpha(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &inflow_base,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Tet4_NS_3D_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  std::string sol_name = Name_Generator(time_info->get_index());
  cur_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  int nl_counter;
  bool rest_flag = restart_init_assembly_flag;

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    nsolver_ptr->GenAlpha_VMS_solve( renew_flag, time_info->get_time(),
        time_info->get_step(), inflow_base, pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, infnbc_part, 
        ebc_part, bc_mat, elementv, elements, quad_v, quad_s, lassem_ptr,
        gassem_ptr, lsolver_ptr, cur_velo, cur_disp, conv_flag, nl_counter );

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index());

    if( time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());
    }

    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }
  delete pre_disp; delete cur_disp; delete pre_velo; delete cur_velo;
}

// EOF
