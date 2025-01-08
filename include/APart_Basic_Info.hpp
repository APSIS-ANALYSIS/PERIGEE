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
    // Static method to get the CPU rank
    static int get_cpu_rank(const std::string &fbasename, const int &in_rank)
    {
      const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

      hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      const int cpu_rank = h5r->read_intScalar("Part_Info", "cpu_rank");

      H5Fclose(file_id);
      return cpu_rank;
    }

    // Static method to get the CPU size
    static int get_cpu_size(const std::string &fbasename, const int &in_rank)
    {
      const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

      hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      const int cpu_size = h5r->read_intScalar("Part_Info", "cpu_size");

      H5Fclose(file_id);
      return cpu_size;
    }

  private:
    // Private constructor to prevent instantiation
    APart_Basic_Info() = delete;
};

#endif
