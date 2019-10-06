#ifndef MESH_NURBS_1PATCH_3D_TH_HPP
#define MESH_NURBS_1PATCH_3D_TH_HPP
// ==================================================================
// Mesh_NURBS_1Patch_3D_TH.hpp
// Single-patch 3D NURBS mesh for Taylor-Hood type element.
//
// Date: June 15 2018
// ==================================================================
#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class Mesh_NURBS_1Patch_3D_TH : public IMesh
{
  public:
    Mesh_NURBS_1Patch_3D_TH( const IMesh * const &mesh_v,
        const IMesh * const &mesh_p );

    virtual ~Mesh_NURBS_1Patch_3D_TH();

  private:
    double hx_max, hy_max, hz_max, hx_min, hy_min, hz_min;
    int s_degree, t_degree, u_degree;
    int nFunc, nElem, nElem_x, nElem_y, nElem_z;
    int nLocBas;
    std::vector<double> hx, hy, hz;
};

#endif
