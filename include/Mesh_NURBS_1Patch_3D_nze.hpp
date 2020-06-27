#ifndef MESH_NURBS_1PATCH_3D_NZE_HPP
#define MESH_NURBS_1PATCH_3D_NZE_HPP
// ============================================================================
// Mesh_NURBS_1Patch_3D_nze.hpp
// This is a global mesh object that contains only nonzero elements.
// 
// The nElem_x, nElem_y, nElem_z, and nElem will be the nonzero-measure
// element's number; the hx, hy, and hz vector will be these nonzero-measure
// elements' measure/length. Using this mesh object, the resulting computing
// mesh will not have zero-measure element.
//
// To use this class, one shall call a regular mesh and IEN constructor, then
// call the IEN_NURBS_1Patch_3D_wPtr( IEN_wzero, mesh_wzero ) to remove the zero
// elements' IEN array. Then call Mesh_NURBS_1Patch_3D_nze( mesh_wzero ) to
// update the mesh object.
//
// Date: March 22 2016
// ============================================================================
#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class Mesh_NURBS_1Patch_3D_nze : public IMesh
{
  public:
    Mesh_NURBS_1Patch_3D_nze( IMesh const * const &inmesh );

    virtual ~Mesh_NURBS_1Patch_3D_nze();

    virtual void print_info() const;

    virtual int get_s_degree() const {return s_degree;}
    virtual int get_t_degree() const {return t_degree;}
    virtual int get_u_degree() const {return u_degree;}

    virtual int get_nFunc_x() const {return nFunc_x;}
    virtual int get_nFunc_y() const {return nFunc_y;}
    virtual int get_nFunc_z() const {return nFunc_z;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem_x() const {return nElem_x;}
    virtual int get_nElem_y() const {return nElem_y;}
    virtual int get_nElem_z() const {return nElem_z;}
   
    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual double get_hx_max() const {return hx_max;}
    virtual double get_hy_max() const {return hy_max;}
    virtual double get_hz_max() const {return hz_max;}

    virtual double get_hx_min() const {return hx_min;}
    virtual double get_hy_min() const {return hy_min;}
    virtual double get_hz_min() const {return hz_min;}

    virtual double get_hx(int ee) const;
    virtual double get_hy(int ee) const;
    virtual double get_hz(int ee) const;

    virtual int get_patch_index() const {return 0;}
    virtual int get_nElem_start() const {return 0;}
    virtual int get_nFunc_start() const {return 0;}

    virtual void get_elem_index( const int &ee, int &ex, 
        int &ey, int &ez ) const;


  private:
    double hx_max, hy_max, hz_max, hx_min, hy_min, hz_min;
    int s_degree, t_degree, u_degree;
    int nFunc_x, nFunc_y, nFunc_z, nFunc;
    int nElem_x, nElem_y, nElem_z, nElem;
    int nLocBas;
    std::vector<double> hx, hy, hz;
};

#endif
