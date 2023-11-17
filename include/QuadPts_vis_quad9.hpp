#ifndef QUADPTS_VIS_QUAD9_HPP
#define QUADPTS_VIS_QUAD9_HPP
// ==================================================================
// QuadPts_vis_quad9.hpp
//
// This is a class of quadrature points for 9 node quad element.
//
// We use [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0].
//        [0.5, 0.0], [1.0, 0.5], [0.5, 1.0], [0.0, 0.5], [0.5, 0.5]
// These are the vertex points for the quad elements.
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_quad9 : public IQuadPts
{
  public:
    QuadPts_vis_quad9() = default;

    virtual ~QuadPts_vis_quad9() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\n===== Visualization Points for Quad9 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================== \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y, so the dim is 2
    virtual int get_dim() const {return 2;}

    // num_pts = num_pts_x x num_pts_y
    virtual int get_num_quadPts() const {return 9;}

    virtual int get_num_quadPts_x() const {return 3;}

    virtual int get_num_quadPts_y() const {return 3;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return 0.5;}

  private:
    const double qp[18] { 0.0, 0.0, 
      1.0, 0.0, 
      1.0, 1.0, 
      0.0, 1.0, 
      0.5, 0.0, 
      1.0, 0.5, 
      0.5, 1.0, 
      0.0, 0.5, 
      0.5, 0.5 };
};

#endif
