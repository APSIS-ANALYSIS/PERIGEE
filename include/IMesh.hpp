#ifndef IMESH_HPP
#define IMESH_HPP
// ==================================================================
// IMesh.hpp
// The Interface for Mesh classes, pure virtual to all different mesh.
// 
// In principle, an unstructured mesh only contains the polynomial 
// degree, the number of nodes, and the number of cells.
//
// The additional get funcitons with _x/y/z is designed for 
// structured mesh.
// 
// The additional get_patch_ is designed for multi-patch geometries.
//
// Date: Sept 23rd 2013.
// ==================================================================
#include "Sys_Tools.hpp"

class IMesh
{
  public:
    IMesh(){};
    
    virtual ~IMesh(){};

    virtual void print_info() const = 0;

    virtual int get_s_degree() const = 0;
    
    virtual int get_t_degree() const = 0;
    
    virtual int get_u_degree() const = 0;
  
    // ------------------------------------------------------------------------ 
    // For unstructured meshes, there is no intrinsic definition of direction
    // we may return the degree in either s, t, or u direction, and the
    // underlying assumption is that get_s_degree(), get_t_degree(), and
    // get_u_degree() return the same value. 
    // ------------------------------------------------------------------------ 
    virtual int get_degree() const { return get_s_degree();}
    
    virtual int get_nFunc() const = 0;
    
    virtual int get_nElem() const = 0;
    
    virtual int get_nLocBas() const = 0;
};

#endif
