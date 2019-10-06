#ifndef APART_BASIC_INFO_HPP
#define APART_BASIC_INFO_HPP
// ==================================================================
// APart_Basic_Info.hpp
// This class stores basic mesh partition info:
// 
// cpu_rank : the id of the cpu
// cpu_size : the total number of the cpu for the simulation
// dual_edge_ncommon : the parameter for mesh partition.
//
// Date: Nov. 8th 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class APart_Basic_Info
{
  public:
    // Constructor: read from h5 file by giving the part file base name
    //              and rank
    APart_Basic_Info( const std::string &fileBaseName, 
        const int &in_cpu_rank );

    virtual ~APart_Basic_Info();

    virtual int get_cpu_rank() const {return cpu_rank;}
    virtual int get_cpu_size() const {return cpu_size;}
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual void print_info() const;
  
  private:
    // Disallow default constructor
    APart_Basic_Info(){};  

    int cpu_rank, cpu_size, dual_edge_ncommon;
};

#endif
