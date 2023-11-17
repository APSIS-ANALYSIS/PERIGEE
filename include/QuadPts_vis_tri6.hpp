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
      for(int ii=0; ii<6; ++ii)
        SYS_T::commPrint("%e, %e, %e, %e \n",
            get_qw(ii), qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
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
