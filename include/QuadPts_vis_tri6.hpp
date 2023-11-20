#ifndef QUADPTS_VIS_TRI6_HPP
#define QUADPTS_VIS_TRI6_HPP
// ==================================================================
// QuadPts_vis_tri6.hpp
//
// This is a class of quadrature points for 6 node triangle element.
//
// We use [0.0 , 0.0], [1.0 , 0.0], [0.0 , 1.0], 
//        [0.5 , 0.0], [0.5 , 0.5], [0.0 , 0.5].
// These are the vertex points for the quadratic triangle elements.
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_tri6 : public IQuadPts
{
  public:
    QuadPts_vis_tri6() = default;

    virtual ~QuadPts_vis_tri6() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\n===== Visualization Points for Tri6 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================= \n");
    }

    // it stores the area coordinate of the quadrature points 
    // in the sequence of r-s-t, so the dim is 3
    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return 6;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return 0.5/6.0;}

  private:
    const double qp[18] { 0.0, 0.0, 1.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.5, 0.0, 0.5,
      0.5, 0.5, 0.0,
      0.0, 0.5, 0.5 };
};

#endif
