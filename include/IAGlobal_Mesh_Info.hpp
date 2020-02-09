#ifndef IAGLOBAL_MESH_INFO_HPP
#define IAGLOBAL_MESH_INFO_HPP
// ==================================================================
// IAGlobal_Mesh_Info.hpp
// This is an interface for a simple class holding Global Mesh Info
// which were writen by preprocessors in GMI group.
//
// Date: Nov. 8th 2013
// ==================================================================
#include "Sys_Tools.hpp"

class IAGlobal_Mesh_Info
{
  public:
    IAGlobal_Mesh_Info(){};
    virtual ~IAGlobal_Mesh_Info(){};

    virtual int get_xdegree() const = 0;
    virtual int get_ydegree() const = 0;
    virtual int get_zdegree() const = 0;

    virtual double get_max_hx() const = 0;
    virtual double get_max_hy() const = 0;
    virtual double get_max_hz() const = 0;

    virtual int get_nElem() const = 0;
    virtual int get_nFunc() const = 0;
    
    virtual int get_nLocBas() const = 0;
    virtual int get_probDim() const = 0;

    virtual int get_elemType() const
    {SYS_T::print_fatal("Error: IAGLobal_Mesh_Info::get_elemType is not implemented. \n"); return -1;}

    virtual void print_info() const = 0;
};

#endif
