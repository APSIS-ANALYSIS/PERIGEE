#include "APart_Basic_Info.hpp"

APart_Basic_Info::APart_Basic_Info( const std::string &fileBaseName, 
    const int &in_cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, in_cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  cpu_rank = h5r->read_intScalar("Part_Info", "cpu_rank");
  cpu_size = h5r->read_intScalar("Part_Info", "cpu_size");
  dual_edge_ncommon = h5r->read_intScalar("Part_Info", "dual_edge_ncommon");

  delete h5r;
  H5Fclose( file_id );
}

APart_Basic_Info::~APart_Basic_Info()
{}

void APart_Basic_Info::print_info() const
{
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Basic Partition Information: \n");
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "cpu_rank: %d \n", cpu_rank);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "cpu_size: %d \n", cpu_size);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "dual_edge_ncommon: %d \n", dual_edge_ncommon);
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

// EOF
