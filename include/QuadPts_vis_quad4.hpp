#ifndef QUADPTS_VIS_QUAD4_HPP
#define QUADPTS_VIS_QUAD4_HPP
// ==================================================================
// QuadPts_vis_quad4.hpp
//
// This is a class of quadrature points for 4 node quad element.
//
// We use [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0].
// These are the vertex points for the quad elements.
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_quad4 : public IQuadPts
{
  public:
    QuadPts_vis_quad4() = default;

    virtual ~QuadPts_vis_quad4() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\n===== Visualization Points for Quad4 ===== \n");
      for(int ii=0; ii<4; ++ii)
        SYS_T::commPrint("%e, %e, %e \n", qw[ii], qp[2*ii], qp[2*ii+1]);
      SYS_T::commPrint("========================================= \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y, so the dim is 2
    virtual int get_dim() const {return 2;}

    // num_pts = num_pts_x x num_pts_y
    virtual int get_num_quadPts() const {return 4;}

    virtual int get_num_quadPts_x() const {return 2;}

    virtual int get_num_quadPts_y() const {return 2;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

  private:
    const double qp[8] { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
    const double qw[4] { 0.5, 0.5, 0.5, 0.5 };
};

#endif
