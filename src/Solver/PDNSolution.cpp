#include "PDNSolution.hpp"

PDNSolution::PDNSolution( const APart_Node * const &pNode )
: dof_num( pNode->get_dof() ),
  nlocal( pNode->get_nlocalnode() * dof_num ),
  nghost( pNode->get_nghostnode() * dof_num )
{
  PetscInt * ifrom = new PetscInt [nghost];

  for(int ii=0; ii<pNode->get_nghostnode(); ++ii)
  {
    for(int jj=0; jj<dof_num; ++jj)
      ifrom[ii*dof_num + jj] = pNode->get_node_ghost(ii) * dof_num + jj;
  }

  VecCreateGhost(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, nghost, ifrom,
      &solution);
   
  delete [] ifrom; ifrom = nullptr;
}

PDNSolution::PDNSolution( const APart_Node * const &pNode,
   const int &input_dof_num )
: dof_num( input_dof_num ),
  nlocal( pNode->get_nlocalnode() * dof_num ),
  nghost( pNode->get_nghostnode() * dof_num )
{
  PetscInt * ifrom = new PetscInt [nghost];

  for(int ii=0; ii<pNode->get_nghostnode(); ++ii)
  {
    for(int jj=0; jj<dof_num; ++jj)
      ifrom[ii*dof_num + jj] = pNode->get_node_ghost(ii) * dof_num + jj;
  }

  VecCreateGhost(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, nghost, ifrom,
      &solution);
   
  delete [] ifrom; ifrom = nullptr;
}

PDNSolution::PDNSolution( const PDNSolution &INPUT )
: dof_num( INPUT.get_dof_num() ),
  nlocal( INPUT.get_nlocal() ),
  nghost( INPUT.get_nghost() )
{
  VecDuplicate(INPUT.solution, &solution);
  VecCopy(INPUT.solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

PDNSolution::PDNSolution( const PDNSolution * const &INPUT_ptr )
: dof_num( INPUT_ptr->get_dof_num() ),
  nlocal( INPUT_ptr->get_nlocal() ),
  nghost( INPUT_ptr->get_nghost() )
{
  VecDuplicate(INPUT_ptr->solution, &solution);
  VecCopy(INPUT_ptr->solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

PDNSolution::~PDNSolution()
{
  VecDestroy(&solution);
}

void PDNSolution::Gen_random()
{
  srand(time(NULL));
  
  PetscScalar * val = new PetscScalar[nlocal];
  PetscInt * idx = new PetscInt[nlocal];
  
  for(int ii=0; ii<nlocal; ++ii) 
  {
    val[ii] = SYS_T::gen_randomD_open(0.0, 1.0);
    idx[ii] = ii;
  }

  VecSetValuesLocal(solution, nlocal, idx, val, INSERT_VALUES);

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  delete [] val; val = nullptr; delete [] idx; idx = nullptr;
}

void PDNSolution::Copy(const PDNSolution &INPUT)
{
  SYS_T::print_fatal_if( dof_num != INPUT.get_dof_num(), "Error: PDNSolution::Copy, dof_num does not match.\n");
  SYS_T::print_fatal_if( nlocal != INPUT.get_nlocal(), "Error: PDNSolution::Copy, nlocal does not match.\n");
  SYS_T::print_fatal_if( nghost != INPUT.get_nghost(), "Error: PDNSolution::Copy, nghost does not match.\n");
  
  VecCopy(INPUT.solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::Copy(const PDNSolution * const &INPUT_ptr)
{
  SYS_T::print_fatal_if( dof_num != INPUT_ptr->get_dof_num(), "Error: PDNSolution::Copy, dof_num does not match.\n");
  SYS_T::print_fatal_if( nlocal != INPUT_ptr->get_nlocal(), "Error: PDNSolution::Copy, nlocal does not match.\n");
  SYS_T::print_fatal_if( nghost != INPUT_ptr->get_nghost(), "Error: PDNSolution::Copy, nghost does not match.\n");
  
  VecCopy(INPUT_ptr->solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::GhostUpdate()
{
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

double PDNSolution::Norm_1() const
{
  double norm;
  VecNorm(solution, NORM_1, &norm);
  return norm;
}

double PDNSolution::Norm_2() const
{
  double norm;
  VecNorm(solution, NORM_2, &norm);
  return norm;
}

double PDNSolution::Norm_inf() const
{
  double norm;
  VecNorm(solution, NORM_INFINITY, &norm);
  return norm;
}

void PDNSolution::PlusAX(const PDNSolution &x, const double &a)
{
  VecAXPY(solution, a, x.solution);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::PlusAX(const PDNSolution * const &x_ptr, const double &a)
{
  VecAXPY(solution, a, x_ptr->solution);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::ScaleValue(const double &val)
{
  VecScale(solution, val);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::GetLocalArray( double * const &local_array ) const
{
  Vec lsol;
  double * array;
  VecGhostGetLocalForm(solution, &lsol);
  VecGetArray(lsol, &array);
  for( int ii=0; ii<(nlocal + nghost); ++ii )
    local_array[ii] = array[ii];
  VecRestoreArray(lsol, &array);
  VecGhostRestoreLocalForm(solution, &lsol);
}

void PDNSolution::GetLocalArray( std::vector<double> &local_array ) const
{
  local_array.resize(nlocal+nghost);
  Vec lsol;
  double * array;
  VecGhostGetLocalForm(solution, &lsol);
  VecGetArray(lsol, &array);
  for( int ii=0; ii<(nlocal + nghost); ++ii )
    local_array[ii] = array[ii];
  VecRestoreArray(lsol, &array);
  VecGhostRestoreLocalForm(solution, &lsol);
}

void PDNSolution::Assembly_GhostUpdate()
{
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::PrintWithGhost() const
{
  Vec lsol;
  double * array;
  VecGhostGetLocalForm(solution, &lsol);
  VecGetArray(lsol, &array);
  for(int i=0; i<nlocal+nghost; ++i)
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%d %g \n",
        i, (double)PetscRealPart(array[i]));
  VecRestoreArray(lsol, &array);
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  VecGhostRestoreLocalForm(solution,&lsol);
}

void PDNSolution::PrintNoGhost() const
{
  VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
}

void PDNSolution::WriteBinary(const char * const &file_name) const
{
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERBINARY);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerBinarySkipInfo(viewer);
  PetscViewerFileSetName(viewer, file_name);
  VecView(solution, viewer);
  PetscViewerDestroy(&viewer);
}

void PDNSolution::ReadBinary(const char * const &file_name) const
{
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name,
      FILE_MODE_READ, &viewer);
  VecLoad(solution, viewer);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
  PetscViewerDestroy(&viewer);
}

// EOF
