#include "PDNSolution_heatEqn.hpp"

PDNSolution_heatEqn::PDNSolution_heatEqn( const class APart_Node * const &pNode,
    const class IALocal_BC * const &locbc, int type )
  : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_ZeroTemp(locbc);
      SYS_T::commPrint("===> Initial solution: Zero temperature for heat equation. \n");
      break;
    case 1:
      Init_OneTemp(locbc);
      SYS_T::commPrint("===> Initial solution: Interior 1.0, boundary 0.0, for heat equation. \n");
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_heatEqn: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1); 
  }
}

PDNSolution_heatEqn::~PDNSolution_heatEqn()
{}

void PDNSolution_heatEqn::Init_ZeroTemp( const class IALocal_BC * const &LBC )
{
  VecSet(solution, 0.0);
  GhostUpdate();
}

void PDNSolution_heatEqn::Init_OneTemp( const class IALocal_BC * const &LBC )
{
  VecSet(solution, 1.0);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  int num = LBC->get_Num_LD(0);
  int * index = new int [num];
  double * value = new double [num];

  for(int ii=0; ii<num; ++ii)
  {
    index[ii] = LBC->get_LDN(0, ii);
    value[ii] = 0.0;
  }
  VecSetValues(solution, num, index, value, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  delete [] index; delete [] value;
}

// EOF
