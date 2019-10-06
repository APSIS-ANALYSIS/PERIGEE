#ifndef MESH_NURBS_1PATCH_3D_HPP
#define MESH_NURBS_1PATCH_3D_HPP

#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class Mesh_NURBS_1Patch_3D : public IMesh
{
  public:
    Mesh_NURBS_1Patch_3D( int in_s_degree, int in_t_degree, int in_u_degree,
        double input_hx_max, double input_hy_max, double input_hz_max,
        double input_hx_min, double input_hy_min, double input_hz_min,
        const std::vector<double> &in_sknot, const std::vector<double> &in_tknot,
        const std::vector<double> &in_uknot );

    // ------------------------------------------------------------------------
    // This is a constructor that does not rely on the knot vector. The user may
    // directly enforce the number of elements, the number of polynomial degree.
    // The number of functions are determined consequently by nFunc = nElem +
    // degree. The purpose of this constructor is to ease test code. The
    // generaetion of a test mesh does not need to rely on the knot vectors.
    // ------------------------------------------------------------------------
    Mesh_NURBS_1Patch_3D( const int &in_s_degree, const int &in_t_degree, 
        const int &in_u_degree,
        const double &input_hx_max, const double &input_hy_max, const double &input_hz_max,
        const double &input_hx_min, const double &input_hy_min, const double &input_hz_min,
        const int &inelemx, const int &inelemy, const int &inelemz );

    virtual ~Mesh_NURBS_1Patch_3D();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return s_degree;}
    virtual int get_t_degree() const {return t_degree;}
    virtual int get_u_degree() const {return u_degree;}

    virtual s_int get_nFunc_x() const {return nFunc_x;}
    virtual s_int get_nFunc_y() const {return nFunc_y;}
    virtual s_int get_nFunc_z() const {return nFunc_z;}
    virtual s_int get_nFunc() const {return nFunc;}

    virtual s_int get_nElem_x() const {return nElem_x;}
    virtual s_int get_nElem_y() const {return nElem_y;}
    virtual s_int get_nElem_z() const {return nElem_z;}
    virtual s_int get_nElem() const {return nElem;}

    virtual s_int get_nElem_x_nz() const {return nElem_x_nz;} 
    virtual s_int get_nElem_y_nz() const {return nElem_y_nz;} 
    virtual s_int get_nElem_z_nz() const {return nElem_z_nz;} 
    virtual s_int get_nElem_nz() const {return nElem_nz;}

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

    virtual int get_patch_index() const {return 0;}
    virtual s_int get_nElem_start() const {return 0;}
    virtual s_int get_nFunc_start() const {return 0;}

    virtual void get_elem_index( const s_int &ee, 
        s_int &ex, s_int &ey, s_int &ez) const;

  private:
    double hx_max, hy_max, hz_max;
    double hx_min, hy_min, hz_min;
    int s_degree, t_degree, u_degree;
    s_int nFunc_x, nFunc_y, nFunc_z, nFunc;
    s_int nElem_x, nElem_y, nElem_z, nElem;
    int nLocBas;
    std::vector<double> hx, hy, hz;

    s_int nElem_x_nz, nElem_y_nz, nElem_z_nz, nElem_nz;
};

#endif
