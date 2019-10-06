#ifndef AGLOBAL_MESH_INFO_MULTIPATCH_NURBS_3D_HPP
#define AGLOBAL_MESH_INFO_MULTIPATCH_NURBS_3D_HPP
// ============================================================================
// AGlobal_Mesh_Info_multiPatch_NURBS_3D.hpp
// This is an instantiation of global mesh info class for multi-patch NURBS 3D
// geometries.
//
// Date: Sept 14 2015
// ============================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "HDF5_PartReader.hpp"

class AGlobal_Mesh_Info_multiPatch_NURBS_3D : public IAGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info_multiPatch_NURBS_3D( const HDF5_PartReader * const &h5reader );

    virtual ~AGlobal_Mesh_Info_multiPatch_NURBS_3D();

    virtual int get_xdegree() const {return xdegree;}
    virtual int get_ydegree() const {return ydegree;}
    virtual int get_zdegree() const {return zdegree;}

    virtual double get_max_hx() const {return hx_max;}
    virtual double get_max_hy() const {return hy_max;}
    virtual double get_max_hz() const {return hz_max;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_probDim() const {return probDim;}
    virtual int get_elemType() const {return elemType;}

    virtual void print() const;

  private:
    int xdegree, ydegree, zdegree;
    double hx_max, hy_max, hz_max, hx_min, hy_min, hz_min;

    int nElem, nFunc;

    int nLocBas, probDim, elemType;
};

#endif
