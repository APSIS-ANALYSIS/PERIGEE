#include "Seg_Sol_Tools.hpp"

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


// EOF
