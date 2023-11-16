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
    QuadPts_vis_hex8();
    
    virtual ~QuadPts_vis_hex8() = default;

    virtual void print_info() const;

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
    double qp[24];
    double qw[8];
};

#endif
