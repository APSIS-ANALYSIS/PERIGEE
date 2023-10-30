#ifndef QUADPTS_VIS_HEX27_HPP
#define QUADPTS_VIS_HEX27_HPP
// ==================================================================
// QuadPts_vis_hex27.hpp
//
// This is a class that stores the visualization sampling points in
// a reference hexahedron.
//
// We use four points at [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
//                       [0.5, 0.0, 0.0], [1.0, 0.5, 0.0], [0.5, 1.0, 0.0], [0.0, 0.5, 0.0], [0.5, 0.0, 1.0], [1.0, 0.5, 1.0], [0.5, 1.0, 1.0], [0.0, 0.5, 1.0],
//                       [0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [1.0, 1.0, 0.5], [0.0, 1.0, 0.5], [0.0, 0.5, 0.5], [1.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 1.0, 0.5],
//                       [0.5, 0.5, 0.0], [0.5, 0.5, 1.0], [0.5, 0.5, 0.5] 
// They are the vertex points of the hexahedron.
// 
// Note: We store them like what we did in QuadPts_Gauss_Hex class, so the dim = 3.       
//
// Date Created: Oct. 30 2023
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_hex27 : public IQuadPts
{
  public:
    QuadPts_vis_hex27();
    virtual ~QuadPts_vis_hex27();

    virtual void print_info() const;

    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return 3;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const 
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    double qp[81];
    double qw[27];
};


#endif
