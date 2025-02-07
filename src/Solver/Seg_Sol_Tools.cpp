#include "Seg_Sol_Tools.hpp"

void SEG_SOL_T::UpdateU( const double &val,
    const PDNSolution * const &usol,
    PDNSolution * const &sol )
{
  // Make sure that the dof of the usol is 3 and of sol is 7
  SYS_T::print_fatal_if(usol->get_dof_num() != 3,
      "Error: SEG_SOL_T::UpdateU the usol vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::UpdateU the sol vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(usol->get_nlocal() / 3 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::UpdateU the usol and sol does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(usol->get_nghost() / 3 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::UpdateU the usol and sol does not match in the number of ghost nodes. \n");

  const int nlocal = usol->get_nlocal() / 3;

  Vec lu, ls;
  double * array_u, * array_s;
  
  VecGhostGetLocalForm(usol->solution, &lu);
  VecGhostGetLocalForm(sol->solution, &ls);

  VecGetArray(lu, &array_u);
  VecGetArray(ls, &array_s);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = ii * 7;
    const int ii3 = ii * 3;

    array_s[ii7]   += val * array_u[ii3];
    array_s[ii7+1] += val * array_u[ii3+1];
    array_s[ii7+2] += val * array_u[ii3+2];
  }

  VecRestoreArray(lu, &array_u);
  VecRestoreArray(ls, &array_s);
  VecGhostRestoreLocalForm(usol->solution, &lu);
  VecGhostRestoreLocalForm(sol->solution, &ls);

  sol->GhostUpdate();
}


void SEG_SOL_T::UpdateV( const double &dt, const double &gamma,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const PDNSolution * const &sol,
    PDNSolution * const &dot_sol )
{
  // Make sure the solution vectors' format is correct.
  SYS_T::print_fatal_if( pre_dot_sol->get_dof_num() !=7,
     "Error: SEG_SOL_T::UpdateV the pre_dot_sol vector dimension is wrong.\n");
  
  SYS_T::print_fatal_if( dot_sol->get_dof_num() !=7,
     "Error: SEG_SOL_T::UpdateV the dot_sol vector dimension is wrong.\n");

  SYS_T::print_fatal_if( pre_sol->get_dof_num() !=7,
     "Error: SEG_SOL_T::UpdateV the pre_sol vector dimension is wrong.\n");

  SYS_T::print_fatal_if( sol->get_dof_num() !=7,
     "Error: SEG_SOL_T::UpdateV the sol vector dimension is wrong.\n");

  const int nlocal = pre_dot_sol -> get_nlocal() / 7;

  Vec dot_n, dot_m, sol_n, sol_m;
  double * adot_n, * adot_m, * asol_n, * asol_m;

  VecGhostGetLocalForm(pre_dot_sol->solution, &dot_n);
  VecGhostGetLocalForm(pre_sol->solution, &sol_n);
  VecGhostGetLocalForm(dot_sol->solution, &dot_m);
  VecGhostGetLocalForm(sol->solution, &sol_m);

  VecGetArray(dot_n, &adot_n); VecGetArray(dot_m, &adot_m);
  VecGetArray(sol_n, &asol_n); VecGetArray(sol_m, &asol_m);

  for(int ii=0; ii<nlocal; ++ii)
  {
    int ii7 = ii * 7;
    adot_m[ii7] = (asol_m[ii7] - asol_n[ii7])/(gamma*dt) + (gamma-1.0)*adot_n[ii7]/gamma;
    adot_m[ii7+1] = (asol_m[ii7+1] - asol_n[ii7+1])/(gamma*dt) + (gamma-1.0)*adot_n[ii7+1]/gamma;
    adot_m[ii7+2] = (asol_m[ii7+2] - asol_n[ii7+2])/(gamma*dt) + (gamma-1.0)*adot_n[ii7+2]/gamma;
  }

  VecRestoreArray(dot_n, &adot_n); VecRestoreArray(dot_m, &adot_m);
  VecRestoreArray(sol_n, &asol_n); VecRestoreArray(sol_m, &asol_m);

  VecGhostRestoreLocalForm(pre_dot_sol->solution, &dot_n);
  VecGhostRestoreLocalForm(pre_sol->solution, &sol_n);
  VecGhostRestoreLocalForm(dot_sol->solution, &dot_m);
  VecGhostRestoreLocalForm(sol->solution, &sol_m);
  
  dot_sol -> GhostUpdate(); 
}


void SEG_SOL_T::PlusAiPV(const double &aa, 
    const double &bb, const double &cc,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 4 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 4,
      "Error: SEG_SOL_T::PlusAiPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 4 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 4 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiPV the solu and step does not match in the number of ghost nodes. \n");

  // Obtain the number of local nodes
  const int nlocal = step->get_nlocal() / 4;

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = ii * 7;
    const int ii4 = ii * 4;
    array_sol[ii7]   += aa * array_step[ii4+1];
    array_sol[ii7+1] += aa * array_step[ii4+2];
    array_sol[ii7+2] += aa * array_step[ii4+3];
    array_sol[ii7+3] += bb * array_step[ii4];
    array_sol[ii7+4] += cc * array_step[ii4+1];
    array_sol[ii7+5] += cc * array_step[ii4+2];
    array_sol[ii7+6] += cc * array_step[ii4+3];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}


void SEG_SOL_T::PlusAiPV(const double &aa, 
    const double &bb, const double &cc,
    const APart_Node * const &pnode,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 4 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 4,
      "Error: SEG_SOL_T::PlusAiPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 4 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 4 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiPV the solu and step does not match in the number of ghost nodes. \n");

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  // Obtain the number of local nodes
  const int nlocal = pnode->get_nlocalnode_solid(); 

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = pnode->get_node_loc_solid(ii) * 7;
    const int ii4 = pnode->get_node_loc_solid(ii) * 4;
    array_sol[ii7]   += aa * array_step[ii4+1];
    array_sol[ii7+1] += aa * array_step[ii4+2];
    array_sol[ii7+2] += aa * array_step[ii4+3];
    array_sol[ii7+3] += bb * array_step[ii4];
    array_sol[ii7+4] += cc * array_step[ii4+1];
    array_sol[ii7+5] += cc * array_step[ii4+2];
    array_sol[ii7+6] += cc * array_step[ii4+3];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}

void SEG_SOL_T::Insert_zero_solid_UV( const APart_Node * const &pnode,
    PDNSolution * const &sol )
{
  // Make sure the dof for sol that 7.
  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::Insert_zero_solid_UPV the solution vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != pnode->get_dof(),
      "Error: SEG_SOL_T::Insert_zero_solid_UPV the solution vector dof number does not match that of APart_Node. \n");

  SYS_T::print_fatal_if(sol->get_nlocal() != pnode->get_nlocalnode() * 7,
      "Error: SEG_SOL_T::Insert_zero_solid_UPV the solution vector dimension does not match that of APart_Node. \n");

  Vec lsol;
  double * array_sol;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGetArray(lsol, &array_sol);

  // Obtain the number of local nodes
  const int nlocal_solid = pnode->get_nlocalnode_solid();

  for(int ii=0; ii<nlocal_solid; ++ii)
  {
    const int ii7 = pnode->get_node_loc_solid(ii) * 7;
    array_sol[ii7]   = 0.0; 
    array_sol[ii7+1] = 0.0;
    array_sol[ii7+2] = 0.0;
     
    array_sol[ii7+4] = 0.0;
    array_sol[ii7+5] = 0.0;
    array_sol[ii7+6] = 0.0;
  }
  
  const int nlocal_fluid = pnode->get_nlocalnode_fluid();
  
  for(int ii=0; ii<nlocal_fluid; ++ii)
  {
    const int ii7 = pnode->get_node_loc_fluid(ii) * 7;
    array_sol[ii7]   = 0.0;
    array_sol[ii7+1] = 0.0;
    array_sol[ii7+2] = 0.0;
  }
  VecRestoreArray(lsol, &array_sol);
  VecGhostRestoreLocalForm(sol->solution, &lsol);

  sol->GhostUpdate(); // update the ghost slots
}

void SEG_SOL_T::PlusAiUPV(const double &aa, 
    const double &bb, const double &cc,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 7 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiUPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 7 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 7 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu and step does not match in the number of ghost nodes. \n");

  // Obtain the number of local nodes
  const int nlocal = step->get_nlocal() / 7;

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = ii * 7;
    array_sol[ii7]   += aa * array_step[ii7+0];
    array_sol[ii7+1] += aa * array_step[ii7+1];
    array_sol[ii7+2] += aa * array_step[ii7+2];
    array_sol[ii7+3] += bb * array_step[ii7+3];
    array_sol[ii7+4] += cc * array_step[ii7+4];
    array_sol[ii7+5] += cc * array_step[ii7+5];
    array_sol[ii7+6] += cc * array_step[ii7+6];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}


void SEG_SOL_T::PlusAiUPV(const double &aa, 
    const double &bb, const double &cc,
    const APart_Node * const &pnode,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 7 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiUPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 7 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 7 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiUPV the solu and step does not match in the number of ghost nodes. \n");

  // Obtain the number of local nodes
  const int nlocal = pnode -> get_nlocalnode_solid(); 

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = pnode->get_node_loc_solid(ii) * 7;
    
    array_sol[ii7]   += aa * array_step[ii7+0];
    array_sol[ii7+1] += aa * array_step[ii7+1];
    array_sol[ii7+2] += aa * array_step[ii7+2];
    array_sol[ii7+3] += bb * array_step[ii7+3];
    array_sol[ii7+4] += cc * array_step[ii7+4];
    array_sol[ii7+5] += cc * array_step[ii7+5];
    array_sol[ii7+6] += cc * array_step[ii7+6];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}


void SEG_SOL_T::PlusAiVPV(const double &aa, 
    const double &bb, const double &cc,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 7 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiVPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 7 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 7 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu and step does not match in the number of ghost nodes. \n");

  // Obtain the number of local nodes
  const int nlocal = step->get_nlocal() / 7;

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = ii * 7;
    array_sol[ii7]   += aa * array_step[ii7+4];
    array_sol[ii7+1] += aa * array_step[ii7+5];
    array_sol[ii7+2] += aa * array_step[ii7+6];
    array_sol[ii7+3] += bb * array_step[ii7+3];
    array_sol[ii7+4] += cc * array_step[ii7+4];
    array_sol[ii7+5] += cc * array_step[ii7+5];
    array_sol[ii7+6] += cc * array_step[ii7+6];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}


void SEG_SOL_T::PlusAiVPV(const double &aa, 
    const double &bb, const double &cc,
    const APart_Node * const &pnode,
    const PDNSolution * const &step,
    PDNSolution * const &sol )
{
  // Make sure the dof for the two vec that 7 and 7.
  SYS_T::print_fatal_if(step->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiVPV the step vector dimension is wrong. \n");

  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu vector dimension is wrong. \n");

  // Make sure the number of local nodes are the same
  SYS_T::print_fatal_if(step->get_nlocal() / 7 != sol->get_nlocal() / 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu and step does not match in the number of local nodes. \n");

  // Make sure the number of ghost nodes are the same
  SYS_T::print_fatal_if(step->get_nghost() / 7 != sol->get_nghost() / 7,
      "Error: SEG_SOL_T::PlusAiVPV the solu and step does not match in the number of ghost nodes. \n");

  // Obtain the number of local nodes
  const int nlocal = pnode->get_nlocalnode_solid();

  Vec lsol, lstep;
  double * array_sol, * array_step;

  VecGhostGetLocalForm(sol->solution, &lsol);
  VecGhostGetLocalForm(step->solution, &lstep);
  VecGetArray(lsol, &array_sol);
  VecGetArray(lstep, &array_step);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii7 = pnode->get_node_loc_solid(ii) * 7;
    array_sol[ii7]   += aa * array_step[ii7+4];
    array_sol[ii7+1] += aa * array_step[ii7+5];
    array_sol[ii7+2] += aa * array_step[ii7+6];
    array_sol[ii7+3] += bb * array_step[ii7+3];
    array_sol[ii7+4] += cc * array_step[ii7+4];
    array_sol[ii7+5] += cc * array_step[ii7+5];
    array_sol[ii7+6] += cc * array_step[ii7+6];
  }

  VecRestoreArray(lsol, &array_sol);
  VecRestoreArray(lstep, &array_step);
  VecGhostRestoreLocalForm(sol->solution, &lsol);
  VecGhostRestoreLocalForm(step->solution, &lstep);

  sol->GhostUpdate(); // update the ghost slots
}

void SEG_SOL_T::Insert_plug_inflow_UPV(const double &val,
    const ALocal_InflowBC * const &infnbc, 
    PDNSolution * const &sol )
{
  // Make sure the dof of the sol vector is 7
  SYS_T::print_fatal_if(sol->get_dof_num() != 7, "Error: SEG_SOL_T::Insert_plug_inflow_UPV the sol vector dimension is wrong. \n");

  SYS_T::print_fatal_if(infnbc->get_num_nbc() != 1, "Error: SEG_SOL_T::Insert_plug_inflow_UPV currently only supports single inlet.\n");

  const int nbc_id = 0;

  const int numnode = infnbc -> get_Num_LD( nbc_id );

  const double dir_x = infnbc -> get_outvec( nbc_id ).x();
  const double dir_y = infnbc -> get_outvec( nbc_id ).y();
  const double dir_z = infnbc -> get_outvec( nbc_id ).z();

  for(int ii=0; ii<numnode; ++ii)
  {
    const int node_index = infnbc -> get_LDN(nbc_id, ii);
    
    VecSetValue(sol->solution, node_index*7+4, dir_x * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*7+5, dir_y * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*7+6, dir_z * val, INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}

void SEG_SOL_T::Extract_solid_U( const APart_Node * const &pnode,
    const PDNSolution * const &sol, PDNSolution * const &output )
{
  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::Extract_solid_U the solution vector dimension is wrong. \n");

  SYS_T::print_fatal_if(output->get_dof_num() != 3,
      "Error: SEG_SOL_T::Extract_solid_U the output vector dimension is wrong. \n");

  Vec lsol;           VecGhostGetLocalForm(sol->solution, &lsol);
  double * array_sol; VecGetArray(lsol, &array_sol);

  Vec lout;           VecGhostGetLocalForm(output->solution, &lout);
  double * array_out; VecGetArray(lout, &array_out);

  const int nlocal_solid = pnode->get_nlocalnode_solid();

  for(int ii=0; ii<nlocal_solid; ++ii)
  {
    const int ii7 = pnode->get_node_loc_solid(ii) * 7;
    const int ii3 = pnode->get_node_loc_solid(ii) * 3;
    array_out[ii3]   = array_sol[ii7];
    array_out[ii3+1] = array_sol[ii7+1];
    array_out[ii3+2] = array_sol[ii7+2];
  }

  VecRestoreArray(lsol, &array_sol);
  VecGhostRestoreLocalForm(sol->solution, &lsol);

  VecRestoreArray(lout, &array_out);
  VecGhostRestoreLocalForm(output->solution, &lout);

  output -> GhostUpdate(); // update the ghost slots
}

void SEG_SOL_T::Extract_solid_P( const APart_Node * const &pnode,
    const PDNSolution * const &sol, PDNSolution * const &output ) 
{ 
  SYS_T::print_fatal_if(sol->get_dof_num() != 7,
      "Error: SEG_SOL_T::Extract_solid_U the solution vector dimension is wrong. \n");

  SYS_T::print_fatal_if(output->get_dof_num() != 1,
      "Error: SEG_SOL_T::Extract_solid_U the output vector dimension is wrong. \n");
  
  Vec lsol;           VecGhostGetLocalForm(sol->solution, &lsol);
  double * array_sol; VecGetArray(lsol, &array_sol);

  Vec lout;           VecGhostGetLocalForm(output->solution, &lout);
  double * array_out; VecGetArray(lout, &array_out);

  const int nlocal_solid = pnode->get_nlocalnode_solid();

  for(int ii=0; ii<nlocal_solid; ++ii)
  { 
    const int jj = pnode->get_node_loc_solid(ii);
    array_out[jj] = array_sol[jj*7 + 3];
  }

  const int nlocal_fluid = pnode->get_nlocalnode_fluid();

  for(int ii=0; ii<nlocal_fluid; ++ii)
  {
    const int jj = pnode->get_node_loc_fluid(ii);
    array_out[jj] = 0.0;
  }

  VecRestoreArray(lsol, &array_sol);
  VecGhostRestoreLocalForm(sol->solution, &lsol);

  VecRestoreArray(lout, &array_out);
  VecGhostRestoreLocalForm(output->solution, &lout);
  
  output -> GhostUpdate(); // update the ghost slots
}

void SEG_SOL_T::CheckUV(const PDNSolution * const &dotU,
    const PDNSolution * const &V )
{
  const int nlgn = dotU->get_nlocal() + dotU->get_nghost();
  const int dofNum = dotU -> get_dof_num();
  double * array_U = new double [nlgn];
  double * array_V = new double [nlgn];

  dotU->GetLocalArray(array_U);
  V->GetLocalArray(array_V);

  for(int ii=0; ii<nlgn/dofNum; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      int index1 = ii * dofNum + jj;
      int index2 = ii * dofNum + jj + 4;
      if( !MATH_T::equals(array_U[index1], array_V[index2], 1.0e-15) )
        PetscPrintf( PETSC_COMM_WORLD,
            "Error: U[%d]=%e and V[%d]=%e, diff is %e, ii=%d, jj=%d \n", 
            index1, array_U[index1], index2, array_V[index2],
            array_U[index1] - array_V[index2], ii, jj );
    }
  }

  delete [] array_U; delete [] array_V;
}

void SEG_SOL_T::Check_dotSol_Sol( const double &dt, const double &gamma,
    const PDNSolution * const &dotSol_n,
    const PDNSolution * const &dotSol_m,
    const PDNSolution * const &Sol_n,
    const PDNSolution * const &Sol_m )
{
  const int nlgn = dotSol_n->get_nlocal() + dotSol_n->get_nghost();
  const int dofNum = dotSol_n -> get_dof_num();

  double * array_ds_n = new double [nlgn];
  double * array_ds_m = new double [nlgn];
  double * array_s_n = new double [nlgn];
  double * array_s_m = new double [nlgn];

  dotSol_n -> GetLocalArray( array_ds_n );
  dotSol_m -> GetLocalArray( array_ds_m );
  Sol_n -> GetLocalArray( array_s_n );
  Sol_m -> GetLocalArray( array_s_m );

  for(int ii=0; ii<nlgn/dofNum; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      const int idx = ii * dofNum + jj;
      double left = array_s_m[idx];
      double righ = array_s_n[idx] + dt * array_ds_n[idx]
        + gamma * dt * (array_ds_m[idx] - array_ds_n[idx]);

      if( !MATH_T::equals(left, righ, 1.0e-15) )
        PetscPrintf( PETSC_COMM_WORLD,
            "Error: SEG_SOL_T::Check_dotSol_Sol error at ii = %d, jj = %d \n",
            ii, jj );
    }
  }

  delete [] array_ds_n; delete [] array_ds_m;
  delete [] array_s_n; delete [] array_s_m;
}

// EOF
