#ifndef QUADPTS_VIS_HEX8_HPP
#define QUADPTS_VIS_HEX8_HPP
// ==================================================================
// QuadPts_vis_hex8.hpp
//
// This is a class that stores the visualization sampling points in
// a reference hexahedron.
//
// We use four points at [0,0,0], [1,0,0], [1,1,0], [0,1,0], [0,0,1], [1,0,1], [1,1,1], [0,1,1].
// They are the vertex points of the hexahedron.
// 
// Note: We store them like what we did in QuadPts_Gauss_Hex class, so the dim = 3.       
//
// Date Created: Oct. 24 2023
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_hex8 final : public IQuadPts
{
  public:
    QuadPts_vis_hex8() = default;
    
    ~QuadPts_vis_hex8() override = default;

    void print_info() const override
    {
      SYS_T::commPrint("\n===== Visualization Points for Hex8 ===== \n");
      IQuadPts::print_info();
      SYS_T::commPrint("========================================= \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y-z, so the dim is 4
    int get_dim() const override {return 3;}

    // num_pts = num_pts_x x num_pts_y x num_pts_z
    int get_num_quadPts() const override {return 8;}

    int get_num_quadPts_x() const override {return 2;}

    int get_num_quadPts_y() const override {return 2;}

    int get_num_quadPts_z() const override {return 2;}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[3*ii+comp];}

    double get_qw(const int &ii) const override {return 0.5;}

  private:
    const double qp[24] {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };
};

#endif
