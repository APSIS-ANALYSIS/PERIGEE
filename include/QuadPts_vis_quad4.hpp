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
      IQuadPts::print_info();
      SYS_T::commPrint("========================================== \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y, so the dim is 2
    virtual int get_dim() const {return 2;}

    // num_pts = num_pts_x x num_pts_y
    virtual int get_num_quadPts() const {return 4;}

    virtual int get_num_quadPts_x() const {return 2;}

    virtual int get_num_quadPts_y() const {return 2;}

    virtual double get_qp(const int &ii, const int &comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(const int &ii) const {return 0.5;}

  private:
    const double qp[8] { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
};

#endif
