#ifndef MESH_NURBS_1PATCH_2D_NZE_HPP
#define MESH_NURBS_1PATCH_2D_NZE_HPP
// ============================================================================
// Mesh_NURBS_1Patch_2D_nze.hpp
// This is a global mesh object that contains only the nonzero elements.
//
// Date: Oct. 5th 2015
// ============================================================================
#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class Mesh_NURBS_1Patch_2D_nze : public IMesh
{
  public:
    Mesh_NURBS_1Patch_2D_nze( IMesh const * const &inmesh );

    virtual ~Mesh_NURBS_1Patch_2D_nze();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return s_degree;}
    virtual int get_t_degree() const {return t_degree;}
    virtual int get_u_degree() const
    {
      std::cerr<<"Error: this is a 2D mesh class, get_u_degree() is invalid! \n";
      return -1;
    }

    virtual s_int get_nFunc_x() const {return nFunc_x;}
    virtual s_int get_nFunc_y() const {return nFunc_y;}
    virtual s_int get_nFunc_z() const
    {
      std::cerr<<"Error: this is a 2D mesh class, get_nFunc_z() is invalid! \n";
      return -1;
    }
    virtual s_int get_nFunc() const {return nFunc;}

    virtual s_int get_nElem_x() const {return nElem_x;}
    virtual s_int get_nElem_y() const {return nElem_y;}
    virtual s_int get_nElem_z() const
    {
      std::cerr<<"Error: this is a 2D mesh class, get_nElem_z() is invalid! \n";
      return -1;
    }
    virtual s_int get_nElem() const {return nElem;}


    virtual int get_nLocBas() const {return nLocBas;}

    virtual double get_hx_max() const {return hx_max;}
    virtual double get_hy_max() const {return hy_max;}
    virtual double get_hz_max() const
    {
      std::cerr<<"Error: this is a 2D mesh class, get_hz_max() is invalid! \n";
      return -1.0;
    }

    virtual double get_hx_min() const {return hx_min;}
    virtual double get_hy_min() const {return hy_min;}
    virtual double get_hz_min() const
    {
      std::cerr<<"Error: this is a 2D mesh class, get_hz_min() is invalid! \n";
      return -1.0;
    }

    virtual double get_hx( s_int ee ) const;
    virtual double get_hy( s_int ee ) const;
    virtual double get_hz( s_int ee ) const;

    virtual int get_patch_index() const {return 0;}
    virtual s_int get_nElem_start() const {return 0;}
    virtual s_int get_nFunc_start() const {return 0;}

    virtual void get_elem_index( const s_int &ee,
        s_int &ex, s_int &ey ) const;

  private:
    double hx_max, hy_max, hx_min, hy_min;
    int s_degree, t_degree;
    int nFunc_x, nFunc_y, nFunc;
    int nElem_x, nElem_y, nElem;
    int nLocBas;
    std::vector<double> hx, hy;
};


#endif
