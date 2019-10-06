#ifndef AGLOBAL_MESH_INFO_1PATCH_NURBS_2D_HPP
#define AGLOBAL_MESH_INFO_1PATCH_NURBS_2D_HPP
// ==================================================================
// AGlobal_Mesh_Info_1Patch_NURBS_2D.hpp
// This is an instantiation of glboal mesh info class which is suited
// for single patch geometry with NURBS 2D geometry
// Date: April 14 2014
// ==================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "HDF5_PartReader.hpp"

class AGlobal_Mesh_Info_1Patch_NURBS_2D : public IAGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info_1Patch_NURBS_2D( const HDF5_PartReader * const &h5reader );
    virtual ~AGlobal_Mesh_Info_1Patch_NURBS_2D();

    virtual int get_xdegree() const {return xdegree;}
    virtual int get_ydegree() const {return ydegree;}
    virtual int get_zdegree() const {return 0;}

    virtual double get_max_hx() const {return hx_max;}
    virtual double get_max_hy() const {return hy_max;}
    virtual double get_max_hz() const {return 0.0;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_probDim() const {return probDim;}
    virtual int get_elemType() const {return elemType;}

    virtual void print() const;
  private:
    int xdegree, ydegree;

    double hx_max, hy_max, hx_min, hy_min;

    int nElem, nElem_x, nElem_y;
    int nFunc, nFunc_x, nFunc_y;

    int nLocBas, probDim, elemType;
};

#endif
