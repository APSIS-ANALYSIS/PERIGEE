#ifndef APART_BASIC_INFO_HPP
#define APART_BASIC_INFO_HPP
// ==================================================================
// APart_Basic_Info.hpp
// This class stores basic mesh partition info:
// 
// cpu_rank : index of the cpu
// cpu_size : total number of cpu's for the simulation
// dual_edge_ncommon : parameter for mesh partition.
//
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
    APart_Basic_Info( const std::string &fbasename, const int &in_rank = 0 );

    virtual ~APart_Basic_Info();

    virtual int get_cpu_rank() const {return cpu_rank;}
    
    virtual int get_cpu_size() const {return cpu_size;}
    
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual void print_info() const;
  
  private:
    // --------------------------------------------------------------
    // Disallow the default constructor
    // --------------------------------------------------------------
    APart_Basic_Info(){};  

    int cpu_rank, cpu_size, dual_edge_ncommon;
};

#endif
