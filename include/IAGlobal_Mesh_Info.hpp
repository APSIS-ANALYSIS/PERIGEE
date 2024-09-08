#ifndef IAGLOBAL_MESH_INFO_HPP
#define IAGLOBAL_MESH_INFO_HPP
// ==================================================================
// IAGlobal_Mesh_Info.hpp
// This is an interface for a simple class holding Global Mesh Info
// which were writen by preprocessors in GMI group.
//
// Date: Nov. 8th 2013
// ==================================================================

class IAGlobal_Mesh_Info
{
  public:
    IAGlobal_Mesh_Info() = default;
    
    virtual ~IAGlobal_Mesh_Info() = default;

    virtual int get_xdegree() const = 0;
    virtual int get_ydegree() const = 0;
    virtual int get_zdegree() const = 0;

    virtual int get_nElem() const = 0;
    virtual int get_nFunc() const = 0;
    
    virtual int get_nLocBas() const = 0;
    virtual int get_probDim() const = 0;

    virtual int get_elemType() const = 0;

    virtual void print_info() const = 0;
};

#endif
