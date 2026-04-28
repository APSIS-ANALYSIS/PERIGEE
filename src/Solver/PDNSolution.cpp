#include "PDNSolution.hpp"
#include <random>

PDNSolution::PDNSolution( const APart_Node * const &pNode )
: PDNSolution( pNode, pNode->get_dof() )
{}

PDNSolution::PDNSolution( const APart_Node * const &pNode, int input_dof_num )
: dof_num( input_dof_num ),
  nlocalnode( pNode->get_nlocalnode() ),
  nghostnode( pNode->get_nghostnode() ),
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
  nlocalnode( INPUT.get_nlocalnode() ),
  nghostnode( INPUT.get_nghostnode() ),
  nlocal( INPUT.get_nlocal() ),
  nghost( INPUT.get_nghost() )
{
  VecDuplicate(INPUT.solution, &solution);
  VecCopy(INPUT.solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

PDNSolution &PDNSolution::operator=( const PDNSolution &INPUT ) noexcept
{
  if( this == &INPUT ) return *this;

  SYS_T::print_fatal_if( dof_num != INPUT.get_dof_num(),
      "Error: PDNSolution::operator=, dof_num does not match.\n" );
  SYS_T::print_fatal_if( nlocal != INPUT.get_nlocal(),
      "Error: PDNSolution::operator=, nlocal does not match.\n" );
  SYS_T::print_fatal_if( nghost != INPUT.get_nghost(),
      "Error: PDNSolution::operator=, nghost does not match.\n" );

  VecCopy(INPUT.solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);

  return *this;
}

PDNSolution::PDNSolution( const PDNSolution * INPUT_ptr )
: dof_num( INPUT_ptr->get_dof_num() ),
  nlocalnode( INPUT_ptr->get_nlocalnode() ),
  nghostnode( INPUT_ptr->get_nghostnode() ),
  nlocal( INPUT_ptr->get_nlocal() ),
  nghost( INPUT_ptr->get_nghost() )
{
  VecDuplicate(INPUT_ptr->solution, &solution);
  VecCopy(INPUT_ptr->solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

PDNSolution::~PDNSolution() noexcept
{
  VecDestroy(&solution);
}

/*
PDNSolution PDNSolution::Gen_random( const APart_Node * const &pNode,
    int input_dof_num )
{
  PDNSolution sol( pNode,
      input_dof_num > 0 ? input_dof_num : pNode->get_dof() );

  thread_local std::mt19937_64 gen( std::random_device{}() );
  std::uniform_real_distribution<double> dis(-1.0, 1.0);

  PetscScalar * val = new PetscScalar[sol.nlocal];
  PetscInt * idx = new PetscInt[sol.nlocal];
  
  for(int ii=0; ii<sol.nlocal; ++ii)
  {
    val[ii] = dis(gen);
    idx[ii] = ii;
  }

  VecSetValuesLocal(sol.solution, sol.nlocal, idx, val, INSERT_VALUES);

  VecAssemblyBegin(sol.solution);
  VecAssemblyEnd(sol.solution);

  sol.GhostUpdate();

  delete [] val; val = nullptr; delete [] idx; idx = nullptr;

  return sol;
}

PDNSolution PDNSolution::Gen_zero( const APart_Node * const &pNode,
    int input_dof_num )
{
  PDNSolution sol( pNode,
      input_dof_num > 0 ? input_dof_num : pNode->get_dof() );

  VecSet(sol.solution, 0.0);

  VecAssemblyBegin(sol.solution);
  VecAssemblyEnd(sol.solution);

  sol.GhostUpdate();

  return sol;
}
*/

void PDNSolution::Copy(const PDNSolution &INPUT)
{
  SYS_T::print_fatal_if( dof_num != INPUT.get_dof_num(), "Error: PDNSolution::Copy, dof_num does not match.\n");
  SYS_T::print_fatal_if( nlocal != INPUT.get_nlocal(), "Error: PDNSolution::Copy, nlocal does not match.\n");
  SYS_T::print_fatal_if( nghost != INPUT.get_nghost(), "Error: PDNSolution::Copy, nghost does not match.\n");
  
  VecCopy(INPUT.solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::Copy(const PDNSolution * INPUT_ptr)
{
  SYS_T::print_fatal_if( dof_num != INPUT_ptr->get_dof_num(), "Error: PDNSolution::Copy, dof_num does not match.\n");
  SYS_T::print_fatal_if( nlocal != INPUT_ptr->get_nlocal(), "Error: PDNSolution::Copy, nlocal does not match.\n");
  SYS_T::print_fatal_if( nghost != INPUT_ptr->get_nghost(), "Error: PDNSolution::Copy, nghost does not match.\n");
  
  VecCopy(INPUT_ptr->solution, solution);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::CopyScale(const PDNSolution &src, double a)
{
  SYS_T::print_fatal_if( !is_layout_equal(*this, src), "Error: PDNSolution::CopyScale, vector layout does not match.\n" );
  SYS_T::print_fatal_if( dof_num != src.get_dof_num(), "Error: PDNSolution::CopyScale, dof_num does not match.\n" );

  VecCopy(src.solution, solution);
  VecScale(solution, a);

  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::CopyScale(const PDNSolution * src_ptr, double a)
{
  CopyScale(*src_ptr, a);
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

void PDNSolution::PlusAX(const PDNSolution &x, double a)
{
  VecAXPY(solution, a, x.solution);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::PlusAX(const PDNSolution * x_ptr, double a)
{
  VecAXPY(solution, a, x_ptr->solution);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::PlusAX(const Vec &x, double a)
{
  VecAXPY(solution, a, x);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::ScaleValue(double val)
{
  VecScale(solution, val);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
}

void PDNSolution::LinearCombination(double a, const PDNSolution &x, double b, const PDNSolution &y)
{
  SYS_T::print_fatal_if( !is_layout_equal(*this, x), "Error: PDNSolution::LinearCombination, vector layout does not match between this and x.\n" );
  SYS_T::print_fatal_if( !is_layout_equal(*this, y), "Error: PDNSolution::LinearCombination, vector layout does not match between this and y.\n" );
  SYS_T::print_fatal_if( !(dof_num == x.get_dof_num() && dof_num == y.get_dof_num()), "Error: PDNSolution::LinearCombination, dof_num does not match.\n" );

  VecCopy(x.solution, solution);
  VecScale(solution, a);
  VecAXPY(solution, b, y.solution);

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

std::vector<double> PDNSolution::GetLocalArray() const
{
  std::vector<double> local_array(nlocal+nghost, 0.0);
  Vec lsol;
  double * array;
  VecGhostGetLocalForm(solution, &lsol);
  VecGetArray(lsol, &array);
  for( int ii=0; ii<(nlocal + nghost); ++ii ) 
    local_array[ii] = array[ii];
  VecRestoreArray(lsol, &array);
  VecGhostRestoreLocalForm(solution, &lsol);

  return local_array;
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

void PDNSolution::WriteBinary(const std::string &file_name) const
{
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERBINARY);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerBinarySkipInfo(viewer);
  PetscViewerFileSetName(viewer, file_name.c_str());
  VecView(solution, viewer);
  PetscViewerDestroy(&viewer);
}

void PDNSolution::ReadBinary(const std::string &file_name)
{
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name.c_str(),
      FILE_MODE_READ, &viewer);
  VecLoad(solution, viewer);
  VecGhostUpdateBegin(solution, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(solution, INSERT_VALUES, SCATTER_FORWARD);
  PetscViewerDestroy(&viewer);
}

bool is_layout_equal( const PDNSolution &left, const PDNSolution &right ) noexcept
{
  return ( left.nlocalnode == right.nlocalnode && left.nghostnode == right.nghostnode );
}

// EOF
