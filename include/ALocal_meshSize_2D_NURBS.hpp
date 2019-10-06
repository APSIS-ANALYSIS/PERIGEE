#ifndef ALOCAL_MESHSIZE_2D_NURBS_HPP
#define ALOCAL_MESHSIZE_2D_NURBS_HPP
// ==================================================================
// ALocal_meshSize_2D_NURBS.hpp
// Local mesh size calss for 2D NURBS mesh
//
// Date:
// April 13 2014
// ==================================================================
#include "HDF5_PartReader.hpp"
#include "IALocal_meshSize.hpp"

class ALocal_meshSize_2D_NURBS : public IALocal_meshSize
{
  public:
    ALocal_meshSize_2D_NURBS(const HDF5_PartReader * const &h5reader);
    virtual ~ALocal_meshSize_2D_NURBS();

    virtual void print() const;

    virtual double get_hx(const int &e) const {return hx[e];}
    virtual double get_hy(const int &e) const {return hy[e];}

    virtual double get_meshsize(const int &e) const {return hx[e] * hy[e];}
  private:
    std::vector<double> hx, hy;
};
#endif
