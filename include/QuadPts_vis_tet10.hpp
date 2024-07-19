#ifndef QUADPTS_VIS_TET10_HPP
#define QUADPTS_VIS_TET10_HPP
// ==================================================================
// QuadPts_vis_tet10.hpp
//
// This is a class that stores the visualization sampling points in
// a reference tetrahedron.
//
// We use [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], 
// [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0], [0, 0, 0.5], 
// [0, 0.5, 0.5], [0.5, 0, 0.5]. 
// They are the vertex points for the quadratic tetrahedron.
//
// Note: We store them in area-coordinates, like what we did in the
//       QuadPts_Gauss_Tet class, so the dim = 4.
//
// Author: Ju Liu
// Date created: March 20 2018
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_tet10 : public IQuadPts
{
  public:
    QuadPts_vis_tet10() = default;

    virtual ~QuadPts_vis_tet10() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\n===== Visualization Points for Tet10 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================= \n");
    }

    virtual int get_dim() const {return 4;}

    virtual int get_num_quadPts() const {return 10;}

    virtual double get_qp(const int &ii, const int &comp) const
    {return qp[4*ii+comp];}

    virtual double get_qw(const int &ii) const {return 0.1/6.0;}

  private:
    const double qp[40] { 0.0, 0.0, 0.0, 1.0,
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      0.5, 0.0, 0.0, 0.5,
      0.5, 0.5, 0.0, 0.0,
      0.0, 0.5, 0.0, 0.5,
      0.0, 0.0, 0.5, 0.5,
      0.0, 0.5, 0.5, 0.0,
      0.5, 0.0, 0.5, 0.0 };
};

#endif
