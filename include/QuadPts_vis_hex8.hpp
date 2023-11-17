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

class QuadPts_vis_hex8 : public IQuadPts
{
  public:
    QuadPts_vis_hex8() = default;
    
    virtual ~QuadPts_vis_hex8() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\n===== Visualization Points for Hex8 ===== \n");
      for(int ii=0; ii<8; ++ii)
        SYS_T::commPrint("%e, %e, %e, %e \n", qw[ii], qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
      SYS_T::commPrint("========================================= \n");
    }

    // it stores the coordinate of the quadrature points 
    // in the sequence of x-y-z, so the dim is 4
    virtual int get_dim() const {return 3;}

    // num_pts = num_pts_x x num_pts_y x num_pts_z
    virtual int get_num_quadPts() const {return 8;}

    virtual int get_num_quadPts_x() const {return 2;}

    virtual int get_num_quadPts_y() const {return 2;}

    virtual int get_num_quadPts_z() const {return 2;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const 
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    const double qp[24] {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

    const double qw[8] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
};

#endif
