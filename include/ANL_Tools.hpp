#ifndef ANL_TOOLS_HPP
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
  int get_int_data(const std::string &fbasename, const int &in_rank, 
      const std::string &partname, const std::string &dataname )
  {
    const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

    hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

    const int val = h5r->read_intScalar(partname.c_str(), dataname.c_str());

    H5Fclose(file_id);
    return val;
  }

  int get_cpu_rank(const std::string &fbasename, const int &in_rank)
  {
    return get_int_data(fbasename, in_rank, "Part_Info", "cpu_rank");
  }

  int get_cpu_size(const std::string &fbasename, const int &in_rank)
  {
    return get_int_data(fbasename, in_rank, "Part_Info", "cpu_size");
  }

  int get_nLocBas(const std::string &fbasename, const int &in_rank)
  {
    return get_int_data(fbasename, in_rank, "Global_Mesh_Info", "nLocBas");
  }

  FE_Typle get_elemType(const std::string &fbasename, const int &in_rank)
  {
    const std::string fName = SYS_T::gen_partfile_name(fbasename, in_rank);

    hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

    auto elemType = FE_T::to_FEType(h5r->read_string("Global_Mesh_Info", "elemType"));

    H5Fclose(file_id);
    return elemType;
  }

} // END OF ANL_T

#endif
