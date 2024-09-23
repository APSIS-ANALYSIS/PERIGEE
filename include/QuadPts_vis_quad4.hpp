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

class QuadPts_vis_quad4 final : public IQuadPts
{
  public:
    QuadPts_vis_quad4() = default;

    ~QuadPts_vis_quad4() override = default;

    void print_info() const override
    {
      SYS_T::commPrint("\n===== Visualization Points for Quad4 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================== \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y, so the dim is 2
    int get_dim() const override {return 2;}

    // num_pts = num_pts_x x num_pts_y
    int get_num_quadPts() const override {return 4;}

    int get_num_quadPts_x() const override {return 2;}

    int get_num_quadPts_y() const override {return 2;}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[2*ii+comp];}

    double get_qw(const int &ii) const override {return 0.5;}

  private:
    const double qp[8] { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
};

#endif
