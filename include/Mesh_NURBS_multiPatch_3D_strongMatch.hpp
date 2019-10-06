#ifndef MESH_NURBS_MULTIPATCH_3D_STRONGMATCH_HPP
#define MESH_NURBS_MULTIPATCH_3D_STRONGMATCH_HPP
// ============================================================================
// Mesh_NURBS_multiPatch_3D_strongMatch.hpp
// 
// This is a mesh handler for multi-patch NURBS geometries. The underlying
// assumption for this mesh object is that the degree for every patch has to be
// identical (hence nLocBas is identical). The maximum/minimum of the mesh size
// will be picked from all the patches.
//
// Date: Aug. 31 2015 
// ============================================================================
#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class Mesh_NURBS_multiPatch_3D_strongMatch : public IMesh
{
  public:
    Mesh_NURBS_multiPatch_3D_strongMatch( 
        const std::vector<IMesh *> &in_mlist, const int &in_num_pat );

    virtual ~Mesh_NURBS_multiPatch_3D_strongMatch();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return s_degree;}
    virtual int get_t_degree() const {return t_degree;}
    virtual int get_u_degree() const {return u_degree;}

    virtual s_int get_nFunc_x() const 
    {std::cerr<<"Error: get_nFunc_x not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual s_int get_nFunc_y() const 
    {std::cerr<<"Error: get_nFunc_y not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual s_int get_nFunc_z() const 
    {std::cerr<<"Error: get_nFunc_z not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    
    virtual s_int get_nFunc() const {return nFunc;}

    virtual s_int get_nElem_x() const
    {std::cerr<<"Error: get_nElem_x not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    virtual s_int get_nElem_y() const
    {std::cerr<<"Error: get_nElem_y not implemented. \n"; exit(EXIT_FAILURE); return 0;}
    virtual s_int get_nElem_z() const
    {std::cerr<<"Error: get_nElem_z not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual s_int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return nLocBas;}
    
    virtual double get_hx_max() const {return hx_max;}
    virtual double get_hy_max() const {return hy_max;}
    virtual double get_hz_max() const {return hz_max;}

    virtual double get_hx_min() const {return hx_min;}
    virtual double get_hy_min() const {return hy_min;}
    virtual double get_hz_min() const {return hz_min;}

    virtual double get_hx( s_int ee ) const;
    virtual double get_hy( s_int ee ) const;
    virtual double get_hz( s_int ee ) const;

    virtual IMesh * get_patch_mesh(const int &pp) const {return mlist[pp];}
    
    virtual int get_num_patch() const {return numPat;}

    // returns the pind index if
    // elem_ptr[pind] <= ee < elem_ptr[pind+1]
    // and 
    // loc_ee = ee - elem_ptr[pind]
    void get_locelem_index( const int &ee, int &pind, int &loc_ee) const;

  private:
    double hx_max, hy_max, hz_max;
    double hx_min, hy_min, hz_min;

    int s_degree, t_degree, u_degree;

    s_int nFunc, nElem;

    int nLocBas;

    // The following are new data structures for the multipatch case
    const unsigned int numPat;
    std::vector<IMesh *> mlist;

    // elem_ptr is an auxilliary array that saves the starting and ending index
    // of element in each patch. It starts with 0. Given the global element
    // index ee, the ee element belongs to patch ii if 
    // elem_ptr[ii] <= ee < elem_ptr[ii+1]
    std::vector<s_int> elem_ptr;
};

#endif
