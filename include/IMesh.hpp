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
    virtual int get_nFunc() const = 0;
    virtual int get_nElem() const = 0;
    virtual int get_nLocBas() const = 0;

    virtual int get_nFunc_x() const
    {SYS_T::print_exit("Error: get_nFunc_x is not implemented. \n"); return 0;}
    
    virtual int get_nFunc_y() const
    {SYS_T::print_exit("Error: get_nFunc_y is not implemented. \n"); return 0;}
    
    virtual int get_nFunc_z() const
    {SYS_T::print_exit("Error: get_nFunc_z is not implemented. \n"); return 0;}
    
    virtual int get_nElem_x() const
    {SYS_T::print_exit("Error: get_nElem_x is not implemented. \n"); return 0;}
    
    virtual int get_nElem_y() const
    {SYS_T::print_exit("Error: get_nElem_y is not implemented. \n"); return 0;}
    
    virtual int get_nElem_z() const
    {SYS_T::print_exit("Error: get_nElem_z is not implemented. \n"); return 0;}
    
    virtual double get_hx_max() const
    {SYS_T::print_exit("Error: get_hx_max is not implemented. \n"); return 0.0;}
    
    virtual double get_hy_max() const
    {SYS_T::print_exit("Error: get_hy_max is not implemented. \n"); return 0.0;}
    
    virtual double get_hz_max() const
    {SYS_T::print_exit("Error: get_hz_max is not implemented. \n"); return 0.0;}

    virtual double get_hx_min() const
    {SYS_T::print_exit("Error: get_hx_min is not implemented. \n"); return 0.0;}
    
    virtual double get_hy_min() const
    {SYS_T::print_exit("Error: get_hy_min is not implemented. \n"); return 0.0;}
    
    virtual double get_hz_min() const
    {SYS_T::print_exit("Error: get_hz_min is not implemented. \n"); return 0.0;}

    virtual double get_hx(int ee) const
    {SYS_T::print_exit("Error: get_hx is not implemented. \n"); return 0.0;}
    
    virtual double get_hy(int ee) const
    {SYS_T::print_exit("Error: get_hy is not implemented. \n"); return 0.0;}
    
    virtual double get_hz(int ee) const
    {SYS_T::print_exit("Error: get_hz is not implemented. \n"); return 0.0;}


    // ------------------------------------------------------------------------
    // Nonzero element numbering info
    // ------------------------------------------------------------------------
    virtual int get_nElem_x_nz() const
    {std::cerr<<"Error: get_nElem_x_nz is not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual int get_nElem_y_nz() const
    {std::cerr<<"Error: get_nElem_y_nz is not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual int get_nElem_z_nz() const
    {std::cerr<<"Error: get_nElem_z_nz is not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual int get_nElem_nz() const
    {std::cerr<<"Error: get_nElem_nz is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    
    // ------------------------------------------------------------------------
    // Information for multipatch
    // ------------------------------------------------------------------------
    virtual int get_patch_index() const
    {std::cerr<<"Error: get_patch_index is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_nElem_start() const
    {std::cerr<<"Error: get_nElem_start is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_nFunc_start() const
    {std::cerr<<"Error: get_nElem_start is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    // ------------------------------------------------------------------------
    // get_patch_mesh belongs to the multipatch mesh handler which returns the
    // invididual mesh's pointer    
    // ------------------------------------------------------------------------
    virtual IMesh * get_patch_mesh(const int &pp) const
    {std::cerr<<"Error: get_patch_mesh is not implmented. \n"; exit(EXIT_FAILURE); return NULL;}

    virtual int get_num_patch() const
    {std::cerr<<"Error: get_num_patch is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    // ------------------------------------------------------------------------
    // get_locelem_index
    // input : global element index ee
    // output: patch index and the element index in this patch
    // ------------------------------------------------------------------------
    virtual void get_locelem_index( const int &ee, int &patch, int &loc_ee) const
    {std::cerr<<"Error: get_locelem_index is not implemented. \n"; exit(EXIT_FAILURE);}
};

#endif
