#include "PETSc_Tools.hpp"

std::string PETSc_T::get_version()
{
  int major, minor, subminor;
  PetscGetVersionNumber(&major,&minor,&subminor, NULL);

  std::string output = "petsc-" + std::to_string( major) + "." + std::to_string( minor)
    + "." + std::to_string( subminor);

  return output;
}

void PETSc_T::MatInfo_Display_local( const Mat &K, 
    const PetscMPIInt &rank )
{
  MatInfo info;
  MatGetInfo(K, MAT_LOCAL, &info);

  PetscMPIInt cpu_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &cpu_rank);

  if(cpu_rank == rank)
  {
    std::cout<<"MatInfo on cpu "<<cpu_rank<<'\n';
    std::cout<<"Block size : "<<info.block_size<<'\n';
    std::cout<<"Nonzero allocated : "<<info.nz_allocated<<'\n';
    std::cout<<"Nonzero used : "<<info.nz_used<<'\n';
    std::cout<<"Nonzero unneeded : "<<info.nz_unneeded<<'\n';
    std::cout<<"Memory : "<<info.memory<<'\n';
    std::cout<<"Number of assemblies : "<<info.assemblies<<'\n';
    std::cout<<"Number of mallocs : "<<info.mallocs<<'\n';
    std::cout<<"Fill ratio given : "<<info.fill_ratio_given<<'\n';
    std::cout<<"Fill raito needed : "<<info.fill_ratio_needed<<'\n';
    std::cout<<"Malloc during factorization : "<<info.factor_mallocs<<'\n';
  }
}

void PETSc_T::MatInfo_Display_global( const Mat &K )
{
  MatInfo info;
  MatGetInfo(K, MAT_GLOBAL_SUM, &info);

  PetscMPIInt cpu_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &cpu_rank);

  if(cpu_rank == 0)
  {
    std::cout<<"MatInfo -- Global Sum \n";
    std::cout<<"Block size : "<<info.block_size<<'\n';
    std::cout<<"Nonzero allocated : "<<info.nz_allocated<<'\n';
    std::cout<<"Nonzero used : "<<info.nz_used<<'\n';
    std::cout<<"Nonzero unneeded : "<<info.nz_unneeded<<'\n';
    std::cout<<"Memory : "<<info.memory<<'\n';
    std::cout<<"Number of assemblies : "<<info.assemblies<<'\n';
    std::cout<<"Number of mallocs : "<<info.mallocs<<'\n';
    std::cout<<"Fill ratio given : "<<info.fill_ratio_given<<'\n';
    std::cout<<"Fill raito needed : "<<info.fill_ratio_needed<<'\n';
    std::cout<<"Malloc during factorization : "<<info.factor_mallocs<<'\n';
  }
}

void PETSc_T::Get_dnz_onz( const Mat &K, 
    std::vector<int> &dnz, std::vector<int> &onz )
{
  PetscInt rstart, rend, cstart, cend;
  MatGetOwnershipRange(K, &rstart, &rend);
  MatGetOwnershipRangeColumn(K, &cstart, &cend);
  const int len = rend - rstart;
  
  dnz.resize(len);
  onz.resize(len);
  
  for(int ii=0; ii<len; ++ii)
  {
    dnz[ii] = 0; onz[ii] = 0;

    PetscInt ncols;
    const PetscInt * cols;
    MatGetRow( K, rstart + ii, &ncols, &cols, NULL );

    for(int jj=0; jj<ncols; ++jj)
    {
      if( cols[jj] >= cstart && cols[jj] < cend ) dnz[ii] += 1;
      else onz[ii] += 1;
    }

    MatRestoreRow(K, rstart + ii, &ncols, &cols, NULL );
  }
}

void PETSc_T::MinusSqrtVec(Vec &v, const double &tol)
{
  PetscInt nn,ii;
  PetscScalar *v1;
  
  VecGetArray(v, &v1);
  VecGetLocalSize(v, &nn);

  for(ii=0; ii<nn; ++ii)
  {
    if(std::abs(v1[ii]) < tol) v1[ii] = 1.0;
    else v1[ii] = 1.0 / PetscSqrtScalar( std::abs(v1[ii]) );
  }

  VecRestoreArray(v, &v1);
}

void PETSc_T::InvAbsVec(Vec &v, const double &tol)
{
  PetscInt nn,ii;
  PetscScalar *v1;
  
  VecGetArray(v, &v1);
  VecGetLocalSize(v, &nn);

  for(ii=0; ii<nn; ++ii)
  {
    if(std::abs(v1[ii]) < tol) v1[ii] = 1.0;
    else v1[ii] = 1.0 / std::abs(v1[ii]);
  }

  VecRestoreArray(v, &v1);
}

void PETSc_T::MatCreateId( Mat &K, const PetscInt &lrow )
{
  MatCreateAIJ(PETSC_COMM_WORLD, lrow, lrow, PETSC_DETERMINE,
     PETSC_DETERMINE, 1, NULL, 0, NULL, &K);

  MatSetOption(K, MAT_SPD, PETSC_TRUE);

  PetscInt mstart, mend;
  MatGetOwnershipRange(K, &mstart, &mend);

  for(PetscInt ii = mstart; ii < mend; ++ii)
    MatSetValue(K, ii, ii, 1.0, INSERT_VALUES);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

double PETSc_T::GetValue( const Vec &a, const int ii )
{
  double val;
  VecGetValues( a, 1, &ii, &val );
  return val;
}

int PETSc_T::GetLocalGhostSize( const Vec &vv )
{
  Vec lsol;
  VecGhostGetLocalForm(vv, &lsol);
  PetscInt NN;
  VecGetSize(lsol, &NN);
  VecGhostRestoreLocalForm(vv, &lsol);

  return static_cast<int>(NN);
}

void PETSc_T::GetLocalArray( const Vec &vv, double * const &vv_array )
{
  const int vv_size = PETSc_T::GetLocalGhostSize( vv );
  
  Vec lsol;
  VecGhostGetLocalForm(vv, &lsol);
  double * array;
  VecGetArray(lsol, &array);
  for( int ii=0; ii<vv_size; ++ii ) vv_array[ii] = array[ii];
  VecRestoreArray(lsol, &array);
  VecGhostRestoreLocalForm(vv, &lsol);
}

std::vector<double> PETSc_T::GetLocalArray( const Vec &vv )
{
  const int vv_size = PETSc_T::GetLocalGhostSize( vv );
  std::vector<double> vv_vector( vv_size, 0.0 );

  Vec lsol;
  VecGhostGetLocalForm(vv, &lsol);
  double * array;
  VecGetArray(lsol, &array);
  for( int ii=0; ii<vv_size; ++ii ) vv_vector[ii] = array[ii];
  VecRestoreArray(lsol, &array);
  VecGhostRestoreLocalForm(vv, &lsol);

  return vv_vector;
}

void PETSc_T::WriteBinary( const Vec &a, const std::string &file_name )
{
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERBINARY);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerBinarySkipInfo(viewer);
  PetscViewerFileSetName(viewer, file_name.c_str());
  VecView(a, viewer);
  PetscViewerDestroy(&viewer);
}

// EOF
