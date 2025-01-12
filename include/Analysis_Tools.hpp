#ifndef ANL_TOLLS_HPP
#define ANL_TOOLS_HPP
// ==================================================================
// ANL_Tools.hpp
//
// cpu_rank : index of the cpu
// cpu_size : total number of cpu's for the simulation
//
// Author: Ju Liu
// Date: Nov. 8th 2013
// ==================================================================
#include "HDF5_Reader.hpp"

namespace ANL_T
{
    int get_cpu_rank(const std::string &fbasename, const int &in_rank)
    {
      const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

      hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      const int cpu_rank = h5r->read_intScalar("Part_Info", "cpu_rank");

      H5Fclose(file_id);
      return cpu_rank;
    }

    int get_cpu_size(const std::string &fbasename, const int &in_rank)
    {
      const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

      hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      const int cpu_size = h5r->read_intScalar("Part_Info", "cpu_size");

      H5Fclose(file_id);
      return cpu_size;
    }

} // END OF ANL_T

#endif
