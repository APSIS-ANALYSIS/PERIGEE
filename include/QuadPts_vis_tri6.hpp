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

class QuadPts_vis_tri6 final : public IQuadPts
{
  public:
    QuadPts_vis_tri6() = default;

    ~QuadPts_vis_tri6() override = default;

    void print_info() const override 
    {
      SYS_T::commPrint("\n===== Visualization Points for Tri6 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================= \n");
    }

    // it stores the area coordinate of the quadrature points 
    // in the sequence of r-s-t, so the dim is 3
    int get_dim() const override {return 3;}

    int get_num_quadPts() const override {return 6;}

    double get_qp(const int &ii, const int &comp) const override 
    {return qp[3*ii+comp];}

    double get_qw(const int &ii) const override {return 0.5/6.0;}

  private:
    const double qp[18] { 0.0, 0.0, 1.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.5, 0.0, 0.5,
      0.5, 0.5, 0.0,
      0.0, 0.5, 0.5 };
};

#endif
