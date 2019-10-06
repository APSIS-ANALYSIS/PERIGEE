#ifndef ALOCAL_MESHSIZE_3D_NURBS_HPP
#define ALOCAL_MESHSIZE_3D_NURBS_HPP
// ==================================================================
// ALocal_meshSize_3D_NURBS.hpp
// Local mesh size class for 3D NURBS mesh.
//
// Date:
// Nov. 11 2013
// ==================================================================
#include <vector>
#include "HDF5_PartReader.hpp"
#include "IALocal_meshSize.hpp"
#include "Vec_Tools.hpp"

class ALocal_meshSize_3D_NURBS : public IALocal_meshSize
{
  public:
    ALocal_meshSize_3D_NURBS(const HDF5_PartReader * const &h5reader);
    virtual ~ALocal_meshSize_3D_NURBS();

    virtual void print() const;

    virtual double get_meshsize(const int &e) const
    {return hx[e] * hy[e] * hz[e];}

    virtual double get_hx(const int &e) const {return hx[e];}
    virtual double get_hy(const int &e) const {return hy[e];}
    virtual double get_hz(const int &e) const {return hz[e];}

  private:
    std::vector<double> hx, hy, hz;
};
#endif
