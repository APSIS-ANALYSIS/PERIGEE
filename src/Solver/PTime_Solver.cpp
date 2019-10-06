#include "PTime_Solver.hpp"

PTime_Solver::PTime_Solver(
    const std::string &input_name, const int &input_record_freq,
    const int &input_renew_tang_freq, const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{
}


PTime_Solver::~PTime_Solver()
{}


std::string PTime_Solver::Name_Generator(const int &counter) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_Solver::Info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "final time: %e \n", final_time);
  PetscPrintf(PETSC_COMM_WORLD, "solution record frequency : %d \n", sol_record_freq);
  PetscPrintf(PETSC_COMM_WORLD, "tangent update frequency over time steps: %d \n", renew_tang_freq);
  PetscPrintf(PETSC_COMM_WORLD, "solution base name: %s \n", pb_name.c_str());
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PTime_Solver::TM_generalized_alpha(
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  cur_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  int nl_counter;
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0)
      renew_flag = true;
    else
      renew_flag = false;

    nsolver_ptr->Gen_alpha_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, bc_part,
        wei_ptr, ele_ptr, lassem_ptr, gassem_ptr, lsolver_ptr, 
        cur_velo, cur_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(),
        time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());
    }
    
    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp;
  delete pre_velo; delete cur_velo;
}


void PTime_Solver::TM_NewtonRaphson(
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  // allocate solution vector in memory
  PDNSolution * curr_disp = new PDNSolution(*init_disp);
  PDNSolution * next_disp = new PDNSolution(*init_disp);
  PDNSolution * step = new PDNSolution(*init_disp);

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  curr_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  int nl_counter;
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0)
      renew_flag = true;
    else
      renew_flag = false;

    nsolver_ptr->NewtonRaphson_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), curr_disp, step,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, bc_part,
        wei_ptr, ele_ptr, lassem_ptr, gassem_ptr, lsolver_ptr, 
        next_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(),
        time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      next_disp->WriteBinary(sol_name.c_str());
    }
    
    curr_disp->Copy(*next_disp);
  }

  delete curr_disp; delete next_disp; delete step;
}


void PTime_Solver::TM_VMS_generalized_alpha(
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  cur_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  int nl_counter;
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0)
      renew_flag = true;
    else
      renew_flag = false;

    nsolver_ptr->Gen_alpha_VMS_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, bc_part,
        wei_ptr, ele_ptr, lassem_ptr, gassem_ptr, lsolver_ptr, 
        cur_velo, cur_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(),
        time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());
    }
    
    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp;
  delete pre_velo; delete cur_velo;
}


void PTime_Solver::TM_VMS_generalized_alpha_noCache_3D(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // save the initial solution
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
    else
      renew_flag = false;

    nsolver_ptr->Gen_alpha_VMS_noCache_3D_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, bc_part,
        wei_ptr, mSize, bs, bt, bu, extractor, lassem_ptr, gassem_ptr, 
        lsolver_ptr, cur_velo, cur_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(),
        time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      cur_disp->WriteBinary(sol_name.c_str());
    }

    pre_disp->Copy(*cur_disp);
    pre_velo->Copy(*cur_velo);
  }

  delete pre_disp; delete cur_disp;
  delete pre_velo; delete cur_velo;
}


void PTime_Solver::TM_NewtonRaphson(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_ElemBC * const &ebc_part,
    const AInt_Weight * const &wei_ptr,
    FEAElement * const &element,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_meshSize * const &mSize,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr ) const
{
  // allocate solution vector in memory
  PDNSolution * curr_disp = new PDNSolution(*init_disp);
  PDNSolution * next_disp = new PDNSolution(*init_disp);
  PDNSolution * step = new PDNSolution(*init_disp);

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  curr_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  bool rest_flag = restart_init_assembly_flag;
  int nl_counter;

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else
      renew_flag = false;

    nsolver_ptr->NewtonRaphson_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), curr_disp, step,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, ebc_part,
        wei_ptr, element, bs, bt, bu, extractor, mSize,
        lassem_ptr, gassem_ptr, lsolver_ptr, next_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(),
        time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      next_disp->WriteBinary(sol_name.c_str());
    }
    
    curr_disp->Copy(*next_disp);
  }

  delete curr_disp; delete next_disp; delete step;
}


void PTime_Solver::TM_NewtonRaphson(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_disp,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_ElemBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    const AInt_Weight * const &wei_ptr,
    FEAElement * const &element,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_meshSize * const &mSize,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr ) const
{
  // allocate solution vector in memory
  PDNSolution * curr_disp = new PDNSolution(*init_disp);
  PDNSolution * next_disp = new PDNSolution(*init_disp);
  PDNSolution * step      = new PDNSolution(*init_disp);

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  curr_disp->WriteBinary(sol_name.c_str());

  bool conv_flag, renew_flag;
  bool rest_flag = restart_init_assembly_flag;
  int nl_counter;

  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag)
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    nsolver_ptr->NewtonRaphson_solve(renew_flag, time_info->get_time(),
        time_info->get_step(), curr_disp, step,
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, ebc_part,
        bc_mat, wei_ptr, element, bs, bt, bu, extractor, mSize,
        lassem_ptr, gassem_ptr, lsolver_ptr, next_disp, conv_flag, nl_counter);

    time_info->TimeIncrement();

    PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index());

    if(time_info->get_index()%sol_record_freq == 0)
    {
      sol_name = Name_Generator( time_info->get_index() );
      next_disp->WriteBinary(sol_name.c_str());
    }
    
    curr_disp->Copy(*next_disp);
  }

  // Clean-up memory
  delete curr_disp; delete next_disp; delete step;
}


void PTime_Solver::TM_GenAlpha(
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &init_acce,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_ElemBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    const AInt_Weight * const &wei_ptr,
    FEAElement * const &element,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_meshSize * const &mSize,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);
  PDNSolution * pre_acce = new PDNSolution(*init_acce);
  PDNSolution * cur_acce = new PDNSolution(*init_acce);

  // save the initial condition
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

    nsolver_ptr->GenAlpha_solve( renew_flag, time_info->get_time(),
        time_info->get_step(), pre_acce, pre_velo, pre_disp, tmga_ptr, 
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, ebc_part, 
        bc_mat, wei_ptr, element, bs, bt, bu, extractor, mSize, lassem_ptr, 
        gassem_ptr, lsolver_ptr, cur_acce, cur_velo,
        cur_disp, conv_flag, nl_counter );

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
    pre_acce->Copy(*cur_acce);
  }

  delete pre_disp; delete cur_disp; 
  delete pre_velo; delete cur_velo; 
  delete pre_acce; delete cur_acce;
}


void PTime_Solver::TM_GenAlpha(
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
    const ALocal_ElemBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    const AInt_Weight * const &wei_ptr,
    FEAElement * const &element,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const IAExtractor * const &extractor,
    const IALocal_meshSize * const &mSize,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // save the initial condition
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

    nsolver_ptr->GenAlpha_solve( renew_flag, time_info->get_time(),
        time_info->get_step(), pre_velo, pre_disp, tmga_ptr, 
        alelem_ptr, lien_ptr, anode_ptr, feanode_ptr, nbc_part, ebc_part, 
        bc_mat, wei_ptr, element, bs, bt, bu, extractor, mSize, lassem_ptr, 
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

  delete pre_disp; delete cur_disp; 
  delete pre_velo; delete cur_velo; 
}


void PTime_Solver::TM_GenAlpha(
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
    PNonlinear_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // save the initial condition
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

    nsolver_ptr->GenAlpha_solve( renew_flag, time_info->get_time(),
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

  delete pre_disp; delete cur_disp; 
  delete pre_velo; delete cur_velo; 
}

// EOF
