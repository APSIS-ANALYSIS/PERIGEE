#ifndef APART_BASIC_INFO_HPP
#define APART_BASIC_INFO_HPP
// ==================================================================
// APart_Basic_Info.hpp
// This class stores basic mesh partition info:
// 
// In the preprocessing stage, there will be the number of CPUs assigned
// for mesh partition; in each file, there will also be a rank ID to
// identify the corresponding CPU of the file. The cpu_size and cpu_rank
// are stored in Part_Info of the hdf5 file. This class will load the two
// from the hdf5 file. 
//
// cpu_rank : index of the cpu
// cpu_size : total number of cpu's for the simulation
//
// Author: Ju Liu
// Date: Nov. 8th 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class APart_Basic_Info
{
  public:
    // --------------------------------------------------------------
    // Constructor: read from the h5 file of the given base name and rank
    //              by default, the complete info is stored on cpu 0.
    //              Therefore, one only needs, and is recommended to,
    //              read from rank 0.
    // --------------------------------------------------------------
    APart_Basic_Info( const std::string &fbasename, const int &in_rank = 0 )
    {
      const std::string fName = SYS_T::gen_partfile_name( fbasename, in_rank );

      hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

      std::unique_ptr<HDF5_Reader> h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      cpu_rank = h5r->read_intScalar("Part_Info", "cpu_rank");
      cpu_size = h5r->read_intScalar("Part_Info", "cpu_size");

      H5Fclose( file_id );
    }

    ~APart_Basic_Info() = default;

    int get_cpu_rank() const {return cpu_rank;}

    int get_cpu_size() const {return cpu_size;}

    void print_info() const
    {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Basic Partition Information: \n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "cpu_rank: %d \n", cpu_rank);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "cpu_size: %d \n", cpu_size);
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }

  private:
    int cpu_rank, cpu_size;

    // Disallow the default constructor
    APart_Basic_Info() = delete; 
};

#endif
